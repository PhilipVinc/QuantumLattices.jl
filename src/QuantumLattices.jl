module QuantumLattices

using Reexport
@reexport using QuantumOptics
import QuantumOptics: basis, SparseOperator, liouvillian
using LightGraphs
using LinearAlgebra, SparseArrays

# 2)
# Then define your abstract types:

abstract type System end
abstract type Problem end

include("Lattices.jl")
include("QuantumHamiltonian.jl")

include("QODefs.jl")
export H_NN_lattice, H_loc_disorder, LossOp_loc, LatticeSumOperator, LatticeHomogeneousState

end # module
