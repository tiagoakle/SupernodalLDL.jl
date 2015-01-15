
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

