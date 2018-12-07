function h_loc_QuantumIsing(; g=0.0)
    return g/2.0*sigmax(SpinBasis(1//2))
end

function h_hop_QuantumIsing(; V=0.0)
    sz = sigmaz(SpinBasis(1//2))
    coeff = V/4.0
    return (sz, sz, coeff)
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

export h_loc_QuantumIsing, h_hop_QuantumIsing, h_loc_tlspumped, h_hop_tightbindingboson
