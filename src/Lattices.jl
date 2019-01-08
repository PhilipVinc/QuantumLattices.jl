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
SquareLattice(dims; PBC=false) = Grid(dims, periodic=PBC)

export SquareLattice
