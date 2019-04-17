
mutable struct GraphLindbladian{G,H,GH,LO} <: AbstractGraphSuperOperator
    graph::G
    hilb::H
    H::GH
    c_ops::LO
end

GraphLindbladian(gr, hilb, Ham) =
    GraphLindbladian(gr, hilb, Ham, Dict{Int,Any}())

function GraphLindbladian(H::AbstractGraphOperator, L=nothing)
    gl = GraphLindbladian(graph(H), basis(H), H)
    !isnothing(L) && add_homogeneous_dissipation(gl, L)
    return gl
end

graph(gl::GraphLindbladian) = gl.graph
basis(gl::GraphLindbladian) = gl.hilb
hamiltonian(gl::GraphLindbladian) = gl.H
jump_operators(gl::GraphLindbladian) = _embed_jump_ops(gl)
liouvillian(gl::GraphLindbladian) = liouvillian(SparseOperator(hamiltonian(gl)),
                                        jump_operators(gl))

function add_homogeneous_dissipation(gl::GraphLindbladian, L)
    !is_homogeneous(gl.hilb) && throw(ErrorException("Can't add homogeneous dissipation to non-homogeneous basis"))

    for v=vertices(graph(gl))
        add_dissipator(gl, L, v)
    end
    return gl
end

function add_dissipator(gl::GraphLindbladian, L, v::Int)
    v ∉ vertices(graph(gl)) && throw(BoundsError("Vertex $v is not a valid site index."))
    basis(gl).bases[v] != basis(L) && throw(ErrorException("For site $v we expect an operator on basis $(basis(gl).bases[v]), but L has $(basis(L))."))
    #TODO Fix multiple loss operators per site
    v ∈ keys(gl.c_ops) && throw(ErrorException("loss op on site $v is already present. Currently we do not yet support multiple loss ops per site."))

    push!(gl.c_ops, v=>L)
    return gl
end

function _embed_jump_ops(gl)
    c_ops = gl.c_ops
    embedded_ops = Vector()

    for v=keys(c_ops)
        op = c_ops[v]
        push!(embedded_ops, embed(basis(gl), v, op))
    end
    return embedded_ops
end

## Integrators
QuantumOptics.steadystate.master(lind::GraphLindbladian, args...) =
    steadystate.master(SparseOperator(hamiltonian(lind)), jump_operators(lind), args...)

QuantumOptics.steadystate.eigenvector(lind::GraphLindbladian, args...) =
    steadystate.eigenvector(SparseOperator(hamiltonian(lind)),
                                            jump_operators(lind), args...)
