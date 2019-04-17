# Given a lindbladian create a local operator
export LocalObservable
"""
    LocalObservable(lindbladian, op)
"""
function LocalObservable(lind::GraphLindbladian, op, norm=1.0)
    Hilb = basis(lind)
    !is_homogeneous(Hilb) && error("The basis should be homogeneous.")
    if isa(op, Function)
        op = op(first(Hilb.bases))
    end

    GraphOperator(graph(lind), Hilb, op/norm)
end
