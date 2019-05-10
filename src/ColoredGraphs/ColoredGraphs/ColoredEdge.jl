using LightGraphs.SimpleGraphs: SimpleEdge, SimpleGraphEdge

import Base: Pair, Tuple, show, ==
import LightGraphs: AbstractEdge, src, dst, reverse
import LightGraphs.SimpleGraphs: SimpleEdge


struct ColoredEdge{T<:Integer} <: AbstractColoredEdge{T}
    src::T
    dst::T
    col::T
end

ColoredEdge(t::Tuple) = ColoredEdge(t[1], t[2], t[3])
ColoredEdge(p::Tuple{Pair,Integer}) = ColoredEdge(p[1].first, p[1].second, p[2])
ColoredEdge(p::Tuple{SimpleEdge,Integer}) = ColoredEdge(src(p[1]), dst(p[1]), p[2])
ColoredEdge(p::SimpleEdge, col) = ColoredEdge(src(p), dst(p), col)
ColoredEdge{T}(p::Tuple{Pair,Integer}) where T<:Integer = ColoredEdge(T(p[1].first), T(p[1].second), T(p[2]))
ColoredEdge{T}(t::Tuple) where T<:Integer = ColoredEdge(T(t[1]), T(t[2]), T(t[3]))
SimpleEdge(ce::AbstractColoredEdge) = SimpleEdge(src(ce), dst(ce))

eltype(e::ET) where ET<:AbstractColoredEdge{T} where T = T

# Accessors
src(e::ColoredEdge) = e.src
dst(e::ColoredEdge) = e.dst
color(e::ColoredEdge) = e.col

# I/O
show(io::IO, e::AbstractColoredEdge) = print(io, "Colored Edge $(e.src) => $(e.dst) = [$(e.col)]")

# Conversions
Pair(e::AbstractColoredEdge) = Pair(src(e), dst(e))
Tuple(e::AbstractColoredEdge) = (src(e), dst(e), color(e))

ColoredEdge{T}(e::AbstractColoredEdge) where T <: Integer = ColoredEdge{T}(T(e.src), T(e.dst),  T(e.col))

# Convenience functions
reverse(e::T) where T<:AbstractColoredEdge = T(dst(e), src(e), color(e))
==(e1::AbstractColoredEdge, e2::AbstractColoredEdge) = (src(e1) == src(e2) && dst(e1) == dst(e2) && color(e1) == color(e2))
