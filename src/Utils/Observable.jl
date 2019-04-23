# Given a lindbladian create a local operator
export LocalObservable
"""
    LocalObservable(lindbladian, op)
"""
function MeanObservable(lind::Union{GraphLindbladian, GraphOperator},
                         op, norm=1.0)
    Hilb = basis(lind)
    !is_homogeneous(Hilb) && error("The basis should be homogeneous.")
    if isa(op, Function)
        op = op(first(Hilb.bases))
    end

    GraphOperator(graph(lind), Hilb, op*norm)
end

function LocalObservable(lind::Union{GraphLindbladian, GraphOperator},
                              op, sites)
    Hilb = basis(lind)
    gop = GraphOperator(graph(lind), Hilb)

    if isa(op, Function) && isa(sites, Int)
        op = op(Hilb.bases[sites])
    elseif isa(op, Function)
    end

    if isa(sites, Int)
        if isa(op, Function)
            op = op(Hilb.bases[sites])
        end
        add_local_operator!(gop, op, sites)
    end
    gop
end
