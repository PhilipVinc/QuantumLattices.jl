@doc raw"""
    H_NN_lattice(lattice, h_loc, hj_1, hj_2, J)

Creates the hamiltonian on `lattice` where h_loc is the local hamiltonian,
`hj_1`and `hj_2` are respectively the jump operator and it's adjointm, and `J`
is a complex-valued proportionality term (defaults to 1.0)
"""
function H_NN_lattice(lattice, h_loc, hj_1, hj_2, J=1.0)
    Hilb = h_loc.basis_l^nv(lattice)

    herm_jumps = ishermitian(hj_1) && ishermitian(hj_2)

    H = SparseOperator(Hilb)
    for v=vertices(lattice)
        # local contribution
        H += embed(Hilb, v, h_loc)

        # NN Hopping (different if hermitian or not)
        J_v  = J*embed(Hilb, v, hj_1)
        if herm_jumps
            for nn=neighbors(lattice, v)
                H += 0.5*J_v*embed(Hilb, nn, hj_2)
            end
        else
            #TODO need to consider the directionality in the case of different
            #jump ops.
            #J_v_dag = conj(J)*embed(Hilb, v, dagger(hj_1))
            #for nn=neighbors(lattice, v)
            #    H += 0.5*J_v*embed(Hilb, nn, hj_2)
            #    H += 0.5*J_v_dag*embed(Hilb, nn, dagger(hj_2))
            #end
        end
    end
    return H
end

function H_loc_disorder(lattice, h_loc, disorder_val)
    Hilb = h_loc.basis_l^nv(lattice)
    H = SparseOperator(Hilb)
    for v=vertices(lattice)
        H += embed(Hilb, v, disorder_val[v]*h_loc)
    end
    return H
end

@doc raw"""
    LossOp_loc(lattice, op::SparseOperator)

Returns a vector with a local-loss operator acting on each lattice site.
"""
function LossOp_loc(lattice, op::SparseOperator)
    LossOp_loc(lattice, [op])
end

@doc raw"""
    LossOp_loc(lattice, ops::Vector{SparseOperator})

Returns a vector with a local-loss operator acting on each lattice site.
"""
function LossOp_loc(lattice, ops::Vector)
    Hilb = ops[1].basis_l^nv(lattice)
    jump_ops = Vector{typeof(embed(Hilb, 1, ops[1]))}()
    for v=vertices(lattice)[1:end]
        for op=ops
            push!(jump_ops, embed(Hilb, v, op))
        end
    end
    return jump_ops
end

@doc raw"""
    LatticeSumOperator(lattice, op::Operator)

Given an operator `op`such as n, returns N=∑ᵢnᵢ for every lattice site.
"""
function LatticeSumOperator(lattice, op)
    Hilb = op.basis_l^nv(lattice)
    opTot = embed(Hilb, 1, op)
    for v=vertices(lattice)[2:end]
        opTot += embed(Hilb, v, op)
    end
    return opTot
end

@doc raw"""
    LatticeHomogeneousState(lattice, ψ_loc)

returns ψ_loc ⊗ ψ_loc ⊗ ... ⊗ ψ_loc for every lattice site.
"""
function LatticeHomogeneousState(lattice, ψ_loc)
    ψ = ψ_loc
    for v=vertices(lattice)[2:end]
        ψ = ψ ⊗ ψ_loc
    end
    return ψ
end
