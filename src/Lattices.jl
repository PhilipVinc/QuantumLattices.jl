export SquareLattice

@doc raw"""
    SquareLattice(dims; PBC=false)

Creates a graph for a Grid-based lattice with dims=[Nx, Ny, Nz...]
dimensions.

Examples:
SquareLattice([10], PBC=true)

creates a PBC chain of 10 sites.

SquareLattice([10, 5], PBC=false)

creates a lattice with 10 sites along x and 5 sites along y, for a total of
50 sites.
"""
SquareLattice(dims; PBC=false) = HyperCube(dims, PBC)

"""
    all_isomorph(graph)

Wraps `LightGraphs.Experimental.all_isomorph` and collects the output in a
format that is compatible with QuantumLattices. Generates a list of all
isomorphisms in the lattice.
"""
function all_isomorph(graph)
    isomorphs = LightGraphs.Experimental.all_isomorph(graph, graph)

    iso_table = Vector{Vector{Int}}()
    for iso=isomorphs
        map_table = zeros(Int, nv(graph))
        for (src,dst)=iso
            map_table[src] = dst
        end
        push!(iso_table, map_table)
    end
    iso_table
end
