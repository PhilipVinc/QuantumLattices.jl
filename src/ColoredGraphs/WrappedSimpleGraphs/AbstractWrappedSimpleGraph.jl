using LightGraphs.SimpleGraphs: SimpleEdge, AbstractSimpleGraph
import LightGraphs.SimpleGraphs: fadj, ne, badj

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

fadj(g::AbstractWrappedSimpleGraph) = fadj(bare(g))
fadj(g::AbstractWrappedSimpleGraph, v::Integer) = fadj(bare(g), v)
badj(g::AbstractWrappedSimpleGraph) = badj(bare(g))
badj(g::AbstractWrappedSimpleGraph, v::Integer) = badj(bare(g), v)

ne(g::AbstractWrappedSimpleGraph) = ne(bare(g))
edgetype(g::AbstractWrappedSimpleGraph) = edgetype(bare(g))
#is_directed(g::AbstractWrappedSimpleGraph) = is_directed(bare(g))
# TODO Fix this ugly hack
is_directed(g::Type{<:AbstractWrappedSimpleGraph}) = is_directed(fieldtype(g, :bare_graph))
has_edge(g::AbstractWrappedSimpleGraph, args...) = has_edge(bare(g), args...)
eltype(g::AbstractWrappedSimpleGraph) = eltype(bare(g))
