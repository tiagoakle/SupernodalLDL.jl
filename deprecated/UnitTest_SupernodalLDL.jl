#Implements tests for the SupernodalLDL code
import Base: SparseMatrixCSC
using Base.Test
include("./SupernodalLDL.jl")

#############################################################
# Tests for IndexUnion function
#############################################################

#Set up a CSR matrix with disjoint non zero sets two nonzeros
#per column and the last column full 
A = spzeros(10,5)
i = 0
for j=1:4
		for k=1:2
		 	i+=1
			A[i,j] = 1
		end
end
for j=1:10
		A[j,5] = 1
end

#Call IndexUnion on the first 4 columns, it should return indices 1,...,8
marks = zeros(Int,10);
(nzr,nnzr) = IndexUnion(A,1,4,marks)
@test nzr[1:nnzr] == [ i for i=1:8]
@test marks == zeros(A.m)

#Call again on the last 3 columns
(nzr,nnzr) = IndexUnion(A,3,3,marks)
@test nzr[1:nnzr] == [i for i=1:10]
@test marks == zeros(A.m)

#Call with a supernode that exceeds the matrix bounds
@test_throws(BoundsError,IndexUnion(A,4,3,marks))

#Call with a supernode that starts after the matrix bounds
@test_throws(BoundsError,IndexUnion(A,6,3,marks))

#Call with an empty supernode
(nzr,nnzr) = IndexUnion(A,3,0,marks)
@test(nnzr==0)

#Call on an empty matrix 
(nzr,nnzr) = IndexUnion(spzeros(10,5),3,3,marks)
@test(nnzr==0)
