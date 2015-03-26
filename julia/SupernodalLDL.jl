import Base: SparseMatrixCSC
#Structure for the supernodal cholesky factor 

typealias LDNumber Union(FloatingPoint, Complex)
#Coordinates: The position of an entry in the full version of the sparse matrix 

#Structures to store a supernodal LDL factor of a Symmetric Matrix
#Each supernode j is a group of s_j contiguous columns with indices c_j,...,c_j+s_j 
type LDFactor{T<:LDNumber}

    n::Int            #Size of the matrix
    ns::Int           #Number of supernodes 
    
    s::Array{Int}     #First ColumnCoordinate in the supernode (ns+1)
    v::Array{T}       #Contains all the nonzero entries.
                      #Packed block by block as follows | diag_1 | rows_1 | diag_2 | rows_2 ... 

    sIxV::Array{Int}   #Index at which each supernode's data starts in array v (ns)
    rCo::Array{Int}    #RowCoordinate of the nonzero rows of every supernode |r11,r12,...|r21,r22,...| where r21 = 5 would imply that the first nonzero offdiagonal row of 
                       #supernode 2 is row 5
    sIxRco::Array{Int} #Index of the first entry of rCo corresponding to each supernode size (ns+1)
    cCoS::Array{Int}   #One entry per column, directs to the supernode [1,1,1,2,2,3,4,4,4] etc

    D::Array{T}      
end

#Given the range of columns that forms a supernode 
#and the csr matrix, this function computes the union 
#of the indices of the nonzero rows of all the columns 
#in the supernode.
# Uses the work vector marks of size n, must be initialized to zeros, will be set to zeros
# at the end
#Returns a tuple (nzr::Array(Int),nnzr::Int) where the array nzr is of size n.
#nzr contains the indices of the non zero rows in the entries nzr[1],...,nzr[nnzr]
#A is a matrix in CSR format 
#sStart is the first column of the supernode 
#width is the number of columns of the supernode
#marks is a work vector of size n set to zero, returns zeros
function IndexUnion(A::SparseMatrixCSC, sStart::Int, width::Int, marks::Array{Int}) 	
	nzr = Array(Int,A.m) #Array of indices
	#Validate input bounds
	if(A.m < sStart+width-1)
			throw(BoundsError("Supernode extends beyond the end of the matrix"))
	end
	nnzr = 0
	for c = sStart:(sStart+width-1) #Iterate over each column of the supernode
		for rp = A.colptr[c]:A.colptr[c+1]-1 #Iterate over the indices of each nonzero row
			if(marks[A.rowval[rp]]==0) #We found a new nonzero
					marks[A.rowval[rp]] = 1
					nnzr += 1
					nzr[nnzr] = A.rowval[rp]
			end
		end
	end

	#Clear all the marks and sort the indices
	for j=1:nnzr
			marks[j] = 0 
	end
	#Sort the section of nzr in place and leave the rest of the vector as is
	sort!(nzr,1,nnzr,QuickSort,Base.Order.ForwardOrdering())
	return (nzr,nnzr)
end

#Based on Jack Poulson's class notes code https://github.com/poulson/ATiNLA/blob/master/Multifrontal.ipynb
type Supernode
    start::Int         # first index of this supernode
    size::Int          # number of indices in this supernode
    struct::Array{Int} # factored structure
    
    children::Array{Supernode}
    has_parent::Bool
    parent::Supernode
   	
	#g: Pruned graph of the supernode
	#indices: 
	#statrt: first index
	#parent: Parent in the elimination tree
    function Supernode(g,indices,start=1,parent=nothing)
        this = new()
        this.has_parent = (parent != nothing)
        if parent != nothing
            this.parent = parent
        end
       	
		#TODO Get rid of this 
        struct_set = Set{Int}()
        
        # Check if this diagonal block is small enough to treat as dense
        size = length(indices)
        if size <= leafSize
            this.start = start
            this.size = size
            # Push in the original structure
            for s=1:this.size
                for jSub in g.adjlist[s]
                    if jSub > this.size
                        push!(struct_set,jSub+(start-1))
                    end
                end
            end
            this.children = Array(Supernode,0)
        else
            ind_copy = copy(indices)
            
            gL, gR, p, p_inv = bisect(g)
            sizeL = length(gL.vertices)
            sizeR = length(gR.vertices)
            this.start = start + sizeL + sizeR
            this.size = size - (sizeL+sizeR)
            
            indL = sub(indices,1:sizeL)
            indR = sub(indices,sizeL+1:sizeL+sizeR)
            indS = sub(indices,sizeL+sizeR+1:size)
            for k=1:this.size
                indS[k] = ind_copy[p[sizeL+sizeR+k]]
            end

            # Push in the original structure
            for k=1:this.size
                s = p[sizeL+sizeR+k]
                for jSub in g.adjlist[s]
                    if jSub > size
                        push!(struct_set,jSub+(start-1))
                    end
                end
            end
            
            this.children = Array(Supernode,2)
            for k=1:sizeL
                indL[k] = ind_copy[p[k]]
            end
            for k=1:sizeR
                indR[k] = ind_copy[p[k+sizeL]]
            end
            this.children[1] = Supernode(gL,indL,start,      this)
            this.children[2] = Supernode(gR,indR,start+sizeL,this)
        end

        # Perform the symbolic factorization now that the children are formed
        # (recall that the original structure has already been inserted)
        for c in this.children
            for s in c.struct
                if s >= start + size
                    push!(struct_set,s)
                end
            end
        end
        # Flatten the set into a vector
        this.struct = Array(Int,length(struct_set))
        k = 1
        for j in struct_set
            this.struct[k] = j
            k += 1
        end
        sort!(this.struct)
        
        this
    end
end
