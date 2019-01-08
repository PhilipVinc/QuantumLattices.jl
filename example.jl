# Imports
using QuantumOptics, QuantumLattices, Statistics, LinearAlgebra

# In this example we look at the fermionized driven-dissipative Bose-Hubabrd model.
# Fermionized in the sense that U->∞ and therefore the local hilbert space has {0,1} photons,

# Parameters
N_sites = 5
J = 20.0
F = 1.0
γ = 1.0
Δ=40.0

# Time evolution Parameters
tlist = collect(range(0.0, stop=100.0, length=1001))
num_traj_base=30

# Generate the lattice
lattice = SquareLattice([N_sites], PBC=true)

# Generate the locally-acting hopping elements (J*a/adag)
h_hop1, h_hop2, J = h_hop_tightbindingboson(J=J)
# Local Hamiltonian for the chosen parameters
hloc    = h_loc_tlspumped(F=F, Δ=Δ)

# Global hamiltonian for the lattice. Beware the - sign because of conventions.
H0 = H_NN_lattice(lattice, hloc, h_hop1, h_hop2, -J)

# Generate a list of loss operator (sqrt(gamma)*a) acting on each site of the lattice
c_ops = LossOp_loc(lattice, sqrt(γ)*destroy(FockBasis(1)))
# Initial Configuration
ψ0 = LatticeHomogeneousState(lattice, fockstate(FockBasis(1), 0))

# Compute the Total number of excitations operator. It's the observable we look at.
Ntot  = LatticeSumOperator(lattice, number(FockBasis(1)))


# time evolution of the ordered system
t, sol = timeevolution.mcwf(tlist, ψ0, H0, c_ops, fout=fout=(t, psi)-> expect(Ntot, psi)/norm(psi)^2)

# we now consider random disorder on Δ with amplitude Δ_disorder_σ
Δ_disorder_σ = 0.5
# Generate the local hamiltonian for Δ=1.0 (n)
hΔ_loc = h_loc_tlspumped(F=0, Δ=1.0)
# Generate the local hamiltonian for F=1.0 (a+adag)
hF_loc = h_loc_tlspumped(F=1.0, Δ=0)

# Generate the amplitudes of the disorder for every site
disorder = Δ_disorder_σ*randn(N_sites)
# Generate the disorder component of the hamiltonian
H_disorder = H_loc_disorder(lattice, hΔ_loc, disorder)
# time evolution of the disordered system
t, sol = timeevolution.mcwf(tlist, ψ0, H0+H_disorder, c_ops, fout=fout=(t, psi)-> expect(Ntot, psi)/norm(psi)^2)
