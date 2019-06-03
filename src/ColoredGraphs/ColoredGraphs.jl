module ColoredGraphs

using LightGraphs

import Base:
    eltype, show, ==, Pair, Tuple, copy, length, issubset, reverse, zero, in, iterate

import LightGraphs:
    _NI, AbstractGraph, AbstractEdge, AbstractEdgeIter,
    src, dst, edgetype, nv, ne, vertices, edges, is_directed,
    has_vertex, has_edge, inneighbors, outneighbors,

    indegree, outdegree, degree, has_self_loops, num_self_loops, insorted

export ColoredGraph, ColoredDiGraph, ColoredEdge

export HyperCube, translational_symm_table

abstract type AbstractColoredEdge{T} <: AbstractEdge{T} end
abstract type AbstractColoredGraph{T} <: AbstractGraph{T} end

include("ColoredGraphs/ColoredEdge.jl")
include("ColoredGraphs/ColoredEdgeIter.jl")
include("ColoredGraphs/ColoredGraph.jl")

include("WrappedSimpleGraphs/AbstractWrappedSimpleGraph.jl")
include("WrappedSimpleGraphs/HyperCube.jl")

end # module
