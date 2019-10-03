using LightGraphs.SimpleGraphs: SimpleEdge, AbstractSimpleGraph
using LightGraphs.SimpleGraphs: SimpleGraphs, fadj, ne, badj, add_edge!, rem_edge!
using LightGraphs: ne, edgetype, has_edge, is_directed

"""
    AbstractWrappedSimpleGraph{T} <: AbstractSimpleGraph{T}

A wrapper type that wraps a SimpleGraph, used to store additional info in it's
type. Mainly used to store information about the spatial embedding (in 1D,2D)
and directions (x,y,z...) of the underlying graph.

Must implement the field `bare_graph::SG` that should store the underlying
graph.
"""
abstract type AbstractWrappedSimpleGraph{T} <: AbstractSimpleGraph{T} end

## imple
bare(g::AbstractWrappedSimpleGraph) = g.bare_graph

SimpleGraphs.fadj(g::AbstractWrappedSimpleGraph) = fadj(bare(g))
SimpleGraphs.fadj(g::AbstractWrappedSimpleGraph, v::Integer) = fadj(bare(g), v)
SimpleGraphs.badj(g::AbstractWrappedSimpleGraph) = badj(bare(g))
SimpleGraphs.badj(g::AbstractWrappedSimpleGraph, v::Integer) = badj(bare(g), v)

LightGraphs.ne(g::AbstractWrappedSimpleGraph) = ne(bare(g))
LightGraphs.edgetype(g::AbstractWrappedSimpleGraph) = edgetype(bare(g))
#is_directed(g::AbstractWrappedSimpleGraph) = is_directed(bare(g))
# TODO Fix this ugly hack
LightGraphs.is_directed(g::Type{<:AbstractWrappedSimpleGraph}) = is_directed(fieldtype(g, :bare_graph))
LightGraphs.has_edge(g::AbstractWrappedSimpleGraph, args...) = has_edge(bare(g), args...)
Base.eltype(g::AbstractWrappedSimpleGraph)   = eltype(bare(g))

SimpleGraphs.add_edge!(g::AbstractWrappedSimpleGraph, e::AbstractEdge) =
    add_edge!(bare(g), e)

SimpleGraphs.rem_edge!(g::AbstractWrappedSimpleGraph, e::AbstractEdge) =
    rem_edge!(bare(g), e)
