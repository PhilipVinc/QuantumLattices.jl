export h_loc_XYZ, h_hop_XYZ, xyz_ham, xyz_lind
@doc raw"""
    h_loc_QuantumIsing(spin=1//2;g=0.0)

Creates the local hamiltonian of a Quantum Ising model defined for `spin` spins
(defaults to $S=1/2$ ).

The Hamiltonian is defined as
```
    h_{loc} = \frac{g}{2}\sigma_x
```
"""
h_loc_XYZ(S=1//2) = SparseOperator(SpinBasis(S))


h_hop_XYZ(S=1//2; Jx=0.0, Jy=0.0, Jz=0.0) = [(sigmax(SpinBasis(S)),sigmax(SpinBasis(S)),Jx),
                                             (sigmay(SpinBasis(S)),sigmay(SpinBasis(S)),Jy),
                                             (sigmaz(SpinBasis(S)),sigmaz(SpinBasis(S)),Jz)]

function xyz_ham(gr::Graph, S=1//2; Jx=0.0, Jy=0.0, Jz=0.0)
    go = GraphOperator(gr, SpinBasis(S), h_loc_XYZ(S))
    for hop=h_hop_XYZ(S, Jx=Jx, Jy=Jy, Jz=Jz)
        add_hop_operator!(go, hop)
    end
    go
end
xyz_lind(gr::AbstractGraph, S=1//2; Jx=0.0, Jy=0.0, Jz=0.0, γ=1.0) = GraphLindbladian(xyz_ham(gr, S, Jx=Jx, Jy=Jy, Jz=Jz),
                                                                    h_loss_QuantumIsing(S,γ=γ))
