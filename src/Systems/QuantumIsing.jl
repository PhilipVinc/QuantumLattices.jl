export h_loc_QuantumIsing, h_hop_QuantumIsing, quantum_ising_ham, quantum_ising_lind
@doc raw"""
    h_loc_QuantumIsing(spin=1//2;g=0.0)

Creates the local hamiltonian of a Quantum Ising model defined for `spin` spins
(defaults to $S=1/2$ ).

The Hamiltonian is defined as
```
    h_{loc} = \frac{g}{2}\sigma_x
```
"""
h_loc_QuantumIsing(S=1//2; g=0.0) = g/2.0*sigmax(SpinBasis(S))


h_hop_QuantumIsing(S=1//2; V=0.0) = (sigmaz(SpinBasis(S)),
                                     sigmaz(SpinBasis(S)),
                                     V/4.0 )

h_loss_QuantumIsing(S=1//2; γ=1.0) = sqrt(γ)*sigmam(SpinBasis(S))

quantum_ising_ham(gr::AbstractGraph, S=1//2; g=0.0, V=0.0) = GraphOperator(gr, SpinBasis(S), h_loc_QuantumIsing(S,g=g), h_hop_QuantumIsing(S,V=V))
quantum_ising_lind(gr::AbstractGraph, S=1//2; g=0.0, V=0.0, γ=1.0) = GraphLindbladian(quantum_ising_ham(gr, S, g=g, V=V), h_loss_QuantumIsing(S,γ=γ))
