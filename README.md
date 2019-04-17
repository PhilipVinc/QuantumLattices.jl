# QuantumLattices.jl

A package for working with Open and Closed quantum systems defined on a lattice.

`QuantumLattices` depends on `QuantumOptics` and on `LightGraphs` for most of
it's functionalities.

The package has two functions:
  - Easily generate hamiltonians on a lattice, by specifing the graph, operators
  on vertices (local terms) and operators on edges (hopping terms).
  - Provide functionality to store Operators as Graphs, sometimes usefull for
  particoular computations
