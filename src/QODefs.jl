import SparseArrays: nnz
export data

data(s::SparseMatrixCSC) = s
data(s::SparseVector) = s
data(s::DataOperator) = s.data
data(s::SuperOperator) = s.data

nnz(s::SparseOperator) = nnz(data(s))
is_homogeneous(b::Basis) = true
is_homogeneous(b::CompositeBasis) = all(y->y==b.bases[1],b.bases)
