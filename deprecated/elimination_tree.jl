import Base: SparseMatrixCSC

#How to compute the elimination tree of a sparse matrix A
type SymbolicAnalysis
    n::Int
    colptr::Vector{Int}  #Pointers to the indices of the columns in rowval
    rowval::Vector{Int}  #Row indices of the non zeros
    parent::Vector{Int}  #Pointer to the parent of the node in the 
                         #elimination tree -1 indicates this is a root
end

#Assume the input matrix A is the upper triangular part of the symmetric matrix A in CSC format
#Produces the SymbolicAnalysis structure containing the elimination tree of A, and the nonzero pattern 
#Of L
function DoSymbolic(A::SparseMatrixCSC)
    #Carry out the symbolic analysis for the matrix A and return the 
    #Symbolic Analysis object

    #Work variables
    p_start::Int = 0
    p_end   = 0::Int
    k       = 0::Int
    p       = 0::Int
    p2      = 0::Int
    n       = A.n::Int;
    
    colptr       = Array(Int,n)  #Pointers to the indices of the columns in rowval
    rowcount     = Array(Int,n)  #Non zero count per row of L
    parent       = Array(Int,n)  #Parent vector of the node
    predecessor  = Array(Int,n)  #Earliest known predecessor in the elimination tree

    #Iterate over every column of A
    for k=1:n
        #Initialize as we go 
        parent[k]      = -1   #First time we touch the parent vector at k
        predecessor[k] = -1

        p_start        = A.colptr[k]
        p_end          = A.colptr[k+1]-1

        for j=p_start:p_end #Iterate over each non zero in the column A(1:k,k)

            i=A.rowval[j]   #Non zero index of the jth non zero if column k 
            if(i < k) #As long as this non zero is above the diagonal
                #Traverse the tree up from the leaves and stop at the root, set the 
                #parent of the root to k and the predecessor array for all the nodes we traverse to k 

                #Skip ahead using the predecessor array 
                p2 = i  #p2 keeps the previous node
                p  = predecessor[i]
                
                while p != -1 && p < k #Traverse the tree up to either a root or node k
                     p2 = p
                     predecessor[p2] = k
                     p  = parent[p]
                end 
                # At this point parent[p2] = -1 or parent[p2] = k
                # if -1 then p2 is a root of the subtree, else this subtree has been connected already
                if(p==-1)
                    parent[p2] = k
                    predecessor[p2] = k
                end
            else # i >= k$ the rest of the non zeros of this column are in the lower triangular part
                j = p_end
                if(i>k)
                    error("The matrix is not upper triangular")
                end
            end
        end
    end
    Sym = SymbolicAnalysis(n,A.colptr,A.rowval,parent)
    return Sym
end

