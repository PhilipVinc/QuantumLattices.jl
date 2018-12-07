module QuantumLattices

using QuantumOptics
using LightGraphs
using LinearAlgebra

# 2)
# Then define your abstract types:
abstract type QuantumSystem end
abstract type System end
abstract type Problem end

include("Lattices.jl")
include("HamiltonianDefinitions.jl")
include("QuantumHamiltonian.jl")

export H_NN_lattice, H_loc_disorder, LossOp_loc, LatticeSumOperator, LatticeHomogeneousState

end # module
