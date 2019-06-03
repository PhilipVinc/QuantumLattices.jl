export GraphOperator, add_local_operator!, add_hop_operator!

mutable struct GraphOperator{G} <: AbstractGraphOperator
    graph::G
    hilb
    LocalOperators::Vector{DataOperator}
    EdgeOpList::Dict{Edge, Vector{Tuple{Number, Int, Int}}}

    HopOperatorsList::Vector{DataOperator}
end

GraphOperator(graph, hilb, loc_op=nothing, hop_op=nothing) =
    GraphOperator(graph, hilb^length(vertices(graph)), loc_op, hop_op)

function GraphOperator(graph, hilb::CompositeBasis, loc_op=nothing, hop_op=nothing)
    loc_ops = [SparseOperator(local_hilb) for local_hilb=hilb.bases ]
    hop_ops = Dict{Edge, Vector{Tuple{Real, Int, Int}}}()
    hop_op_list = Vector{DataOperator}()
    go = GraphOperator{typeof(graph)}(graph, hilb, loc_ops, hop_ops, hop_op_list)
    !isnothing(loc_op) && add_local_operator!(go, loc_op)
    !isnothing(hop_op) && add_hop_operator!(go, hop_op)
    return go
end

graph(go::GraphOperator) = go.graph
basis(go::GraphOperator) = go.hilb
data(go::GraphOperator) = data(SparseOperator(go))

function add_local_operator!(go::GraphOperator, loc_ops, sites=vertices(graph(go)))
    if isa(loc_ops, AbstractOperator)
        loc_ops = fill(loc_ops, length(sites))
    end

    for (i,op) = zip(sites, loc_ops)
        go.LocalOperators[i] .+= op
    end
end

add_hop_operator!(go, (hl, hr, coeff)) = add_hop_operator!(go, hl, hr, coeff)
function add_hop_operator!(go::GraphOperator, hl, hr, coeff::Real)
    #imag(coeff) != 0 && !isa(graph(go), DiGraph) && error("If you use complex hopping phase then you should use DiGraphs.")

    herm_jumps = ishermitian(hl) && ishermitian(hr)
    #coeff = herm_jumps ? coeff/2 : coeff

    # Find the id for this operator. If not present, add it.
    l_id = findfirst((x)->x==hl, go.HopOperatorsList)
    if isnothing(l_id)
        push!(go.HopOperatorsList, hl)
        l_id = length(go.HopOperatorsList)
    end
    r_id = findfirst((x)->x==hr, go.HopOperatorsList)
    if isnothing(r_id)
        push!(go.HopOperatorsList, hr)
        r_id = length(go.HopOperatorsList)
    end

    # For every edge store a triplet with the coefficient and the two indices
    # identifying the operator.
    for edge=edges(graph(go))
        edge_op = Tuple{Number, Int, Int}([coeff, l_id, r_id])
        if edge ∉ keys(go.EdgeOpList)
            go.EdgeOpList[edge] = [edge_op]
        else
            push!(go.EdgeOpList[edge], edge_op)
        end
    end
end

function SparseOperator(go::GraphOperator)
    hilb = go.hilb
    op = SparseOperator(hilb)

    # Row, column, value vectors for easier initialization of the operator.
    rows = Vector{Int}()
    cols = Vector{Int}()
    vals = Vector{eltype(first(go.LocalOperators).data)}()
    n_rows, n_cols = size(op.data)

    for site=vertices(graph(go))
        eop = embed(hilb, site, go.LocalOperators[site])
        e_row, e_col, e_val = findnz(eop.data)
        append!(rows, e_row)
        append!(cols, e_col)
        append!(vals, e_val)
        #op .+= eop
    end

    # Cache holding the embedded local operators that form a jump.
    lop_cache = Dict{Tuple{Int,Int},DataOperator}()
    for edge=edges(graph(go)) # for each edge
        if edge ∈ keys(go.EdgeOpList)
            connections = go.EdgeOpList[edge]
            for (coeff, l_id, r_id)=connections
                hl = get!(lop_cache,
                          tuple(edge.src, l_id),
                          embed(hilb, edge.src, go.HopOperatorsList[l_id]))
                hr = get!(lop_cache,
                          tuple(edge.dst, r_id),
                          embed(hilb, edge.dst, go.HopOperatorsList[r_id]))

                # if it's the same operator left and right
                if l_id == r_id
                    eop = coeff* hl*hr
                else # it's different operators
                    eop = coeff*hl*hr +conj(coeff)*hl*hr
                end
                #op .+= eop
                e_row, e_col, e_val = findnz(eop.data)
                append!(rows, e_row)
                append!(cols, e_col)
                append!(vals, e_val)
            end
        end
    end

    op.data = sparse(rows, cols, vals, n_rows, n_cols)

    op
end

# pretty printing
