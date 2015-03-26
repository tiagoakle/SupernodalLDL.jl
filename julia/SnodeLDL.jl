# This module packages Jack Poulson's Supernodal Cholesky code 
# from his Advanced Topics in Numerical Linear Algebra 
# https://github.com/poulson/ATiNLA.
#
# Minor modifications were done for compatibility with Julia 3.6
# Comments added by packagers

using Graphs, Metis

#Transforms a matrix in CSC format to a graph represented by  
#an adjecency structure from Graphs.jl 
function sparseToAdj(M)
    if M.m != M.n
        error("Expected a symmetric matrix")
    end
    n = M.n
    g = simple_adjlist(n,is_directed=false)
    for j=1:n
        for k=M.colptr[j]:M.colptr[j+1]-1
            i = M.rowval[k]
            if i > j
                add_edge!(g,i,j)
            end
        end
    end
    g
end

#Given an array *part* of size n, where
#part[i] in {0,1,2} generate a permutation p 
#such that part[p[i]] is [0,,...,0,1,...,1,2,...,2]
function partToPerm(part)
    sizes   = zeros(Int,3)
    offsets = zeros(Int,3)
    n = length(part)
    for j=1:n
        sizes[part[j]+1] += 1
    end
    offsets[1] = 0
    offsets[2] = sizes[1]
    offsets[3] = sizes[1] + sizes[2]
    
    p = zeros(Int,n)
    for j=1:n
        partVal = part[j]
        p[offsets[partVal+1]+1] = j
        offsets[partVal+1] += 1
    end
    
    p
end

#Given a graph represented as an undirected adjacency list *g* 
#  --possibly with edges leading outside its vertex set-- [1,...,g.vertices]
#1- Partition the vertex set into 3 subsets L,R,S using Metis. 
#2- Generate a permutation p such that p[i] for i in 1,...,|L| is in L,
#   p[i] for i in |L|+1,...,|L|+|R| is in R and the remaining in S.
#   and the inverse permutation p_inv
#3- Generate adjacency lists for the subgraph induced by L and R 
#  	*gL* and *gR*   
#  generate a permutation p and inverse permutation p_inv such that 
#
function bisect(g)
    n = length(g.vertices)
    gPruned = simple_adjlist(n,is_directed=false)
    for s=1:n
        for t in g.adjlist[s]
            if t > s && t <= n
                add_edge!(gPruned,s,t)
            end
        end
    end 
    sizes, part = vertexSep(gPruned)
    p = partToPerm(part)
    p_inv = invperm(p)
    
	sizeL = int(sizes[1])
    sizeR = int(sizes[2])
    
    gL = simple_adjlist(int(sizeL),is_directed=true)
    gR = simple_adjlist(int(sizeR),is_directed=true)
    for s=1:n
        if part[s] == 0
            sMap = int(p_inv[s])
            for t in g.adjlist[s]
                if t <= n
                    tMap = int(p_inv[t])
                    add_edge!(gL,sMap,tMap)
                else
                    add_edge!(gL,sMap,t)
                end
            end
        elseif part[s] == 1
            sMap = int(p_inv[s]-sizeL)
            for t in g.adjlist[s]
                if t <= n
                    tMap = int(p_inv[t]-sizeL)
                    add_edge!(gR,sMap,tMap)
                else
                    add_edge!(gR,sMap,t-sizeL)
                end
            end
        end
    end
    
    gL, gR, p, p_inv
end

type Supernode
    start::Int         # first index of this supernode
    size::Int          # number of indices in this supernode
    struct::Array{Int} # factored structure
    
    children::Array{Supernode}
    has_parent::Bool
    parent::Supernode
    
    function Supernode(g,indices,start=1,parent=nothing)
        this = new()
        this.has_parent = (parent != nothing)
        if parent != nothing
            this.parent = parent
        end
        
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

type Front{F}
    L::Array{F,2}
    BR::Array{F,2}
    
    children::Array{Front{F}}
    child_maps::Array{Any,1}
    has_parent::Bool
    parent::Front{F}
    
    function Front(APerm::SparseMatrixCSC{F},snode::Supernode,parent=nothing)
        this = new()
        this.has_parent = (parent != nothing)
        if this.has_parent
            this.parent = parent
        end
        
        # Initialize this node of the matrix
        structSize = length(snode.struct)
        this.L = zeros(F,snode.size+structSize,snode.size)
        for jSub=1:snode.size
            j = jSub + (snode.start-1)
            for k=APerm.colptr[j]:APerm.colptr[j+1]-1
                i = APerm.rowval[k]
                if i >= snode.start && i < snode.start+snode.size
                    iSub = i - (snode.start-1)
                    this.L[iSub,jSub] = APerm.nzval[k]
                    this.L[jSub,iSub] = conj(APerm.nzval[k])
                elseif i >= snode.start+snode.size
                    structRan = searchsorted(snode.struct,i)
                    if length(structRan) == 1
                        iSub = structRan[1] + snode.size
                        this.L[iSub,jSub] = APerm.nzval[k]
                    else
                        error("Found ",length(structRan)," instances of ($i,$j) in struct(",
                        snode.start,":",snode.start+snode.size-1,")=",snode.struct)
                    end
                end
            end
        end
        this.BR = zeros(F,structSize,structSize)
        
        # Handle the children
        num_children = length(snode.children)
        this.children = Array(Front{F},num_children)
        this.child_maps = Array{Any,1}[]
        for c=1:num_children
            this.children[c] = Front{F}(APerm,snode.children[c],this)
            push!(this.child_maps,Array(Int,length(snode.children[c].struct)))
            for k=1:length(snode.children[c].struct)
                i = snode.children[c].struct[k]
                if i < snode.start+snode.size
                    this.child_maps[c][k] = i - (snode.start-1)
                else
                    loc = searchsorted(snode.struct,i)
                    if length(loc) == 1
                        this.child_maps[c][k] = loc[1] + snode.size
                    else
                        error("expected a single index but found ",loc)
                    end
                end
            end
        end

        this
    end 
end

function Unpack!{F}(front::Front{F},snode::Supernode,L::SparseMatrixCSC{F})
    num_children = length(front.children)
    for c=1:num_children
        Unpack!(front.children[c],snode.children[c],L) 
    end
    
    start = snode.start
    totalSize, nodeSize = size(front.L)
    diagInd = start:start+nodeSize-1
    struct = snode.struct
    structSize = length(struct)

    if nodeSize > 0
        L[diagInd,diagInd] = copy(front.L[1:nodeSize,:])
        if structSize > 0
            L[struct,diagInd] = copy(front.L[nodeSize+1:end,:])
        end
    end
end

function Unpack{F}(front::Front{F},snode::Supernode)
    n = snode.start + snode.size-1
    L = speye(n,n)
    Unpack!(front,snode,L)
    L
end

function Cholesky!{F}(front::Front{F})
    # Recurse on the children
    num_children = length(front.children)
    for c=1:num_children
        Cholesky!(front.children[c])
    end
    
    m,nodeSize = size(front.L)
    structSize = m - nodeSize
    FTL = sub(front.L,1:nodeSize,1:nodeSize)
    FBL = sub(front.L,nodeSize+1:m,1:nodeSize)
    FBR = front.BR
    
    # Perform the extend-adds
    for c=1:num_children
        childStructSize = length(front.child_maps[c]);
        for jChild=1:childStructSize
            jSub = front.child_maps[c][jChild];
            for iChild=jChild:childStructSize
                iSub = front.child_maps[c][iChild];
                value = front.children[c].BR[iChild,jChild];
                if iSub <= nodeSize
                    FTL[iSub,jSub] += value;
                elseif jSub <= nodeSize
                    FBL[iSub-nodeSize,jSub] += value;
                else
                    FBR[iSub-nodeSize,jSub-nodeSize] += value;
                end
            end
        end
        # TODO: Clear front.children[c].BR
	end

    LAPACK.potrf!(string(:L)[1], FTL)
    BLAS.trsm!('R','L','T','N',1.,FTL,FBL)
    BLAS.syrk!('L','N',-1.,FBL,1.,FBR)
end

type RHSNode{F}
    T::Array{F,2}
    B::Array{F,2}
    
    children::Array{RHSNode{F}}
    has_parent::Bool
    parent::RHSNode{F}
    
    function RHSNode(X,snode::Supernode,parent=nothing)
        this = new()
        this.has_parent = (parent != nothing)
        if this.has_parent
            this.parent = parent
        end
        
        n,numRHS = size(X)
        this.T = copy(X[snode.start:snode.start+snode.size-1,:])
        this.B = zeros(F,0,0)
        
        num_children = length(snode.children)
        this.children = Array(RHSNode{F},num_children)
        for c=1:num_children
            this.children[c] = RHSNode{F}(X,snode.children[c],this)
        end
        
        this
    end
end

function Unpack!{F}(XNode::RHSNode{F},snode::Supernode,X::StridedMatrix{F})
    ind = snode.start:snode.start+snode.size-1
    X[ind,:] = copy(XNode.T)
    
    num_children = length(XNode.children)
    for c=1:num_children
        Unpack!(XNode.children[c],snode.children[c],X) 
    end
end

function Unpack{F}(XNode::RHSNode{F},snode::Supernode)
    n = snode.start+snode.size-1
    rootSize, numRHS = size(XNode.T)
    X = zeros(F,n,numRHS)
    Unpack!(XNode,snode,X)
    X
end

function ForwardSolve!{F}(front::Front{F},snode::Supernode,XNode::RHSNode{F})
    # Recurse on the children
    num_children = length(front.children)
    for c=1:num_children
        ForwardSolve!(front.children[c],snode.children[c],XNode.children[c])
    end
    
    # Accumulate the updates from the children
    totalSize, nodeSize = size(front.L)
    nodeSize, numRHS = size(XNode.T)
    structSize = totalSize - nodeSize
    XNode.B = zeros(structSize,numRHS)
    for c=1:num_children
        childStructSize = length(front.child_maps[c]);
        for iChild=1:childStructSize
            iSub = front.child_maps[c][iChild];
            for j=1:numRHS
                value = XNode.children[c].B[iChild,j];
                if iSub <= nodeSize
                    XNode.T[iSub,j] += value;
                else
                    XNode.B[iSub-nodeSize,j] += value;
                end
            end
        end
        # TODO: Clear XNode.children[c].B
    end
    
    LT = sub(front.L,1:nodeSize,1:nodeSize)
    LB = sub(front.L,nodeSize+1:totalSize,1:nodeSize)
    BLAS.trsm!('L','L','N','N',1.,LT,XNode.T)
    BLAS.gemm!('N','N',-1.,LB,XNode.T,1.,XNode.B)
end

function BackwardSolve!{F}(front::Front{F},snode::Supernode,XNode::RHSNode{F})
    totalSize, nodeSize = size(front.L)
    nodeSize, numRHS = size(XNode.T)
    structSize = totalSize - nodeSize
    XNode.B = zeros(structSize,numRHS)
    
    # Pull updates down from the parent
    if front.has_parent
        parentNodeSize = snode.parent.size
        num_siblings = length(front.parent.children)
        
        # Determine which sibling we are
        whichSibling = -1
        for s=1:num_siblings
            if front === front.parent.children[s]
                whichSibling = s
            end
        end
        if whichSibling == -1
            error("This front was not a child of the parent")
        end
        
        for iSub=1:structSize
            iParent = front.parent.child_maps[whichSibling][iSub];
            for j=1:numRHS
                value = XNode.B[iSub,j];
                if iParent <= parentNodeSize
                    XNode.B[iSub,j] += XNode.parent.T[iParent,j]
                else
                    XNode.B[iSub,j] += XNode.parent.B[iParent-parentNodeSize,j]
                end
            end
        end
        # TODO: Clear XNode.parent.B
    end
    
    LT = sub(front.L,1:nodeSize,1:nodeSize)
    LB = sub(front.L,nodeSize+1:totalSize,1:nodeSize)
    BLAS.gemm!('T','N',-1.,LB,XNode.B,1.,XNode.T)
    BLAS.trsm!('L','L','T','N',1.,LT,XNode.T)
    
    # Recurse on the children
    num_children = length(front.children)
    for c=1:num_children
        BackwardSolve!(front.children[c],snode.children[c],XNode.children[c])
    end
end

function Solve{F}(front::Front{F},root::Supernode,p::Array{Int},B::StridedMatrix{F})
    # P A P^T = L L^H and A x = b imply
    #
    #     x = inv(P) (L^H \ (L \ (P b)))
    #
    XNodal = RHSNode{F}(B[p,:],root)
    ForwardSolve!(front,root,XNodal)
    BackwardSolve!(front,root,XNodal)
    XPerm = Unpack(XNodal,root)
    X = XPerm[invperm(p),:]
    X
end
