using Base.Test

include("elimination_tree.jl")
#Tests for the code in the elimination tree 
import Base: spdiagm

#The Cholesky factorization in the Julia code does not return explicit factors in 
#Julia sparse format, therefore we will use a simple inefficient factorization to compute the factors and elimination tree
function inefficient_cholesky(A::SparseMatrixCSC{Float64})
    #Quick checks 
    if(A.m != A.n)
        #Throw exception
        error("Matrix must be square")
    end
    #Threshold for conditioning check
    seps = sqrt(eps(Float64))

    #Matrix for the cholesky factor
    L = tril(A)
    
    #Matrix for the diagonal
    D = Array(Float64,A.n)
   for k=1:n-1
     #Populate the diagonal 
     D[k] = L[k,k]
     d    = D[k]
     if(abs(d)<seps)
         #Throw exception
         error("Matrix has ill conditioned minors")
     end
    
     #Set the diagonal to 1
     L[k,k] = 1;
     #Scale the column
     L[k+1:n,k] = L[k+1:n,k]/d
     #Subtract the outer product
     L[k+1:n,k+1:n] = L[k+1:n,k+1:n] - tril(L[k+1:n,k]*L[k+1:n,k]'*d);
   end
   D[n]   = L[n,n]
   L[n,n] = 1 

   nzix = 0
   colp = 0
   
   #At this point the factorization is complete so we can compute the elimination tree
   T = Array(Int,A.n)
   T[1:n] = -1;
   for k=1:n-1
    colp = L.colptr[k]
    colpn = L.colptr[k+1] 
    if(colpn==colp+1) #ony one non zero, it has to be the diagonal and this is a root
        #Do nothing
    else
       T[k] = L.rowval[colp+1]  #The parent is the first non zero after the diagonal
    end
   end
   return (L,D,T)
end


#Build a sparse upper triangular matrix and do the 
#symbolic analysis. 
srand(0)
n = 40
A = sprand(n,n,0.5)::SparseMatrixCSC{Float64}
A = A'*A+speye(n)

(L,D,T) = inefficient_cholesky(A)
Res = L*spdiagm(D)*L'-A
#First test the above function to have some certainty that our tests will be sensical
@test_approx_eq_eps(maxabs(Res),0,1E-14)

#Use DoSymbolic to generate the tree
Sym = DoSymbolic(triu(A))
@test(isequal(T,Sym.parent))

#Test with an empty matrix 

#Test with a matrix with zero diagonals but that fills in as the factorization proceeds

#Test with a matrix with a zero row
