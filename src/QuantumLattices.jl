module QuantumLattices

using Reexport, Requires
@reexport using QuantumOpticsBase
import QuantumOpticsBase: basis, SparseOperator, liouvillian
import QuantumOpticsBase: DenseOperator, dense
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
export add_local_operator!, add_hop_operator!
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

function __init__()
    @require QuantumOptics="6e0679c1-51ea-5a7c-ac74-d61b76210b0c" begin
        ## Integrators
        QuantumOptics.steadystate.master(lind::GraphLindbladian, args...) =
            QuantumOptics.steadystate.master(SparseOperator(hamiltonian(lind)), jump_operators(lind), args...)

        QuantumOptics.steadystate.eigenvector(lind::GraphLindbladian, args...) =
            QuantumOptics.steadystate.eigenvector(SparseOperator(hamiltonian(lind)),
                                                    jump_operators(lind), args...)
    end

end

end # module
