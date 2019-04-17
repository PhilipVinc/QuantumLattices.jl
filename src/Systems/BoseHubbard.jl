
@doc raw"""
    h_loc_BoseHubbard(n_max; U=0.0, F=0.0, Δ=0.0)

Creates the local hamiltonian of a driven Bose-Hubbard model defined on a
Fock Space with at most `n_max` bosons.

The Hamiltonian is defined as
```
    h_{loc} = -\delta\hat{a}^\dagger\hat{a} + U\hat{a}^\dagger\hat{a}^\dagger\hat{a}\hat{a} + (F^\star\hat{a} + F\hat{a}^\dagger)
```
"""
function h_loc_BoseHubbard(n_max;U=0.0, F=0.0, Δ=0.0)
    Hilb_loc = FockBasis(n_max)
    a_loc  = destroy(Hilb_loc)
    ad_loc = create(Hilb_loc)
    n_loc  = ad_loc*a_loc
    return - Δ*n_loc + U*n_loc*(n_loc-identityoperator(Hilb_loc)) + conj(F)*a_loc + F*ad_loc
end

function h_loc_tlspumped(;F=0.0, Δ=0.0)
    Hilb_loc = FockBasis(1)
    a_loc  = destroy(Hilb_loc)
    ad_loc = create(Hilb_loc)
    n_loc  = ad_loc*a_loc
    return F*(a_loc+ad_loc) - Δ*n_loc
end

function h_hop_tightbindingboson(; J=0.0)
    Hilb_loc = FockBasis(1)
    a_loc  = destroy(Hilb_loc)
    ad_loc = create(Hilb_loc)
    return (a_loc, ad_loc, J)
end

export h_loc_BoseHubbard, h_loc_tlspumped, h_hop_tightbindingboson
