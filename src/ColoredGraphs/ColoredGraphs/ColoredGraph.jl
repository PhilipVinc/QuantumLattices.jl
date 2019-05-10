using LightGraphs.SimpleGraphs: SimpleEdge, AbstractSimpleGraph

mutable struct ColoredGraph{SGT<:AbstractSimpleGraph,T} <: AbstractColoredGraph{T}
    uncolored_graph::SGT
    colmap::Dict{SimpleEdge{T}, T}
end

function ColoredGraph(g::AbstractSimpleGraph)
    T=eltype(g)
    cmap = Dict{SimpleEdge{T}, T}()
    for e=edges(g)
        cmap[e] = 1
    end
    return ColoredGraph(g, cmap)
end

ColoredGraph(n::T=0) where T = ColoredGraph(SimpleGraph(n), Dict{SimpleEdge{T}, T}())
ColoredDiGraph(n::T=0) where T = ColoredGraph(SimpleDiGraph(n), Dict{SimpleEdge{T}, T}())

eltype(x::ColoredGraph{SGT,T}) where {SGT,T} = T

LightGraphs.nv(g::AbstractColoredGraph) = nv(g.uncolored_graph)
LightGraphs.ne(g::AbstractColoredGraph) = ne(g.uncolored_graph)
LightGraphs.vertices(g::AbstractColoredGraph) = vertices(g.uncolored_graph)
LightGraphs.edges(g::AbstractColoredGraph) = ColoredEdgeIter(g)

colmap(g::AbstractColoredGraph) = g.colmap
uncolored(g::AbstractColoredGraph) = g.uncolored_graph

set_color!(g::AbstractColoredGraph, src, dst, col) = begin
    colmap(g)[SimpleEdge(src,dst)] = col
    !is_directed(g) ? colmap(g)[SimpleEdge(dst, src)] = col : nothing;
    return nothing
end
color(g::AbstractColoredGraph, src, dst) = color(g, SimpleEdge(src,dst))
color(g::AbstractColoredGraph, e) = colmap(g)[e]
colored(g::AbstractColoredGraph, e::SimpleEdge) = ColoredEdge(e, color(g, e))
colored(g::AbstractColoredGraph, src, dst) = ColoredEdge(src, dst, color(g, src, dst))

LightGraphs.has_vertex(g::AbstractColoredGraph, v::Integer) = v in vertices(g)

# vertex properties
LightGraphs.inneighbors(g::AbstractColoredGraph, v::Integer) = inneighbors(uncolored(g),v)
LightGraphs.outneighbors(g::AbstractColoredGraph, v::Integer) = outneighbors(uncolored(g),v)

## specific implementations
LightGraphs.is_directed(g::ColoredGraph) = is_directed(uncolored(g))
LightGraphs.has_edge(g::ColoredGraph, s, d) = has_edge(uncolored(g), s, d)

function LightGraphs.has_edge(g::ColoredGraph, s, d, c)
    has_edge(uncolored(g), s, d) || return false
    return color(g, s, d) == c
end
LightGraphs.has_edge(g::ColoredGraph, e::ColoredEdge) = has_edge(g, src(e), dst(e), col(e))

# modifications
LightGraphs.add_edge!(g::AbstractColoredGraph, e::AbstractColoredEdge) =
    add_edge!(g, src(e), dst(e), color(e))
LightGraphs.add_edge!(g::AbstractColoredGraph, src, dst, col) = begin
    !add_edge!(uncolored(g), src, dst) && return false
    set_color!(g, src, dst, col)
    return true
end

LightGraphs.rem_edge!(g::AbstractColoredGraph, e) = rem_edge!(g, src(e), dst(e))
LightGraphs.rem_edge!(g::AbstractColoredGraph, src, dst) = begin
    !rem_edge!(uncolored(g))  && return false
    delete!(colmap(g), SimpleEdge(src,dst))
    !is_directed(g) && delete!(colmap(g), SimpleEdge(dst, src))
    return true
end

LightGraphs.add_vertex!(g::AbstractColoredGraph) = add_vertex!(uncolored(g))
LightGraphs.add_vertices!(g::AbstractColoredGraph, n) = add_vertices!(uncolored(g), n)

LightGraphs.rem_vertex!(g::AbstractColoredGraph, v) = begin
    dsts = outneighbors(g, v)
    !rem_vertex!(uncolored(g)) && return false
    for dst=dsts
        delete!(colmap(g), SimpleEdge(v,dst))
    end
end

# I/O
show(io::IO, g::ColoredGraph) = print(io, "Colored Graph with $(nv(g)) vertices, $(ne(g)) edges")
