using LightGraphs.SimpleGraphs: SimpleEdgeIter, SimpleEdgeIterState, edge_start, edge_next

"""
    ColoredEdgeIter

The function [`edges`](@ref) returns a `SimpleEdgeIter` for `AbstractSimpleGraphs`.
The iterates are in lexicographical order, smallest first. The iterator is valid for
one pass over the edges, and is invalidated by changes to the graph.

# Examples
```jldoctest
julia> using LightGraphs

julia> g = PathGraph(3);

julia> es = edges(g)
SimpleEdgeIter 2

julia> e_it = iterate(es)
(Edge 1 => 2, SimpleEdgeIterState [2, 2])

julia> iterate(es, e_it[2])
(Edge 2 => 3, SimpleEdgeIterState [0, 1])
```
"""
struct ColoredEdgeIter{G} <: AbstractEdgeIter
    g::G
end

eltype(::Type{SimpleEdgeIter{CGT}}) where {CGT<:AbstractColoredGraph} = ColoredEdge{T}
#eltype(::Type{ColoredEdgeIterState{ColoredDiGraph{T}}}) where {T} = ColoredDiGraphEdge{T}

function iterate(eit::ColoredEdgeIter{G}) where {G<:AbstractColoredGraph}
    state = edge_start(uncolored(eit.g))
    return iterate(eit, state)
end

function iterate(eit::ColoredEdgeIter{G}, state::SimpleEdgeIterState{T}) where {T,G<:AbstractColoredGraph{T}}
    state.s == zero(T) && return nothing
    edge, state = edge_next(uncolored(eit.g), state)
    return colored(eit.g, edge), state
end

length(eit::ColoredEdgeIter) = ne(eit.g)

function _isequal(e1::ColoredEdgeIter, e2)
    k = 0
    for e in e2
        has_edge(e1.g, e) || return false
        k += 1
    end
    return k == ne(e1.g)
end

in(e, es::ColoredEdgeIter) = has_edge(es.g, e)

show(io::IO, eit::ColoredEdgeIter) = write(io, "ColoredEdgeIter $(ne(eit.g))")
