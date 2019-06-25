module QuantumLattices

using Reexport
@reexport using QuantumOptics
import QuantumOptics: basis, SparseOperator, liouvillian
#using QuantumOptics : sigmax, sigmay, sigmaz, SpinBasis, FockBasis, destroy, create, identityoperator
#using QuantumOptics : conj, ishermitian, SparseOperator, embed
#
using LinearAlgebra, SparseArrays

# 2)
# Then define your abstract types:

abstract type System end
abstract type Problem end
using LightGraphs
include("ColoredGraphs/ColoredGraphs.jl")
using .ColoredGraphs
export HyperCube, translational_symm_table


include("Lattices.jl")
include("QuantumHamiltonian.jl")

abstract type AbstractGraphOperator end
abstract type AbstractGraphSuperOperator end

export graph, basis, data, is_homogeneous
export add_local_operator!, add_hop_operator!, SparseOperator, DenseOperator
export GraphOperator, GraphLindbladian, hamiltonian, jump_operators, liouvillian
export add_homogeneous_dissipation, add_dissipator

include("QODefs.jl")
include("GraphOperator.jl")
include("GraphSuperOperator.jl")
include("Utils/Observable.jl")

# Prefedined physical systems
include("Systems/QuantumIsing.jl")
include("Systems/BoseHubbard.jl")
include("Systems/XYZ.jl")

export H_NN_lattice, H_loc_disorder, LossOp_loc, LatticeSumOperator, LatticeHomogeneousState

end # module
