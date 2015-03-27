

function ForwardSolveLDL!{F}(front::Front{F},snode::Supernode,XNode::RHSNode{F})
    # Recurse on the children
    num_children = length(front.children)
    for c=1:num_children
        ForwardSolveLDL!(front.children[c],snode.children[c],XNode.children[c])
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
    BLAS.trsm!('L','L','N','U',1.,LT,XNode.T)
    BLAS.gemm!('N','N',-1.,LB,XNode.T,1.,XNode.B)
end

function BackwardSolveLDL!{F}(front::Front{F},snode::Supernode,XNode::RHSNode{F})
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
	#DIVIDE by diagonal
	D  = diag(front.L);
	XNode.T = XNode.T./D
    BLAS.trsm!('L','L','T','N',1.,LT,XNode.T)
    
    # Recurse on the children
    num_children = length(front.children)
    for c=1:num_children
        BackwardSolveLDL!(front.children[c],snode.children[c],XNode.children[c])
    end
end

function SolveLDL{F}(front::Front{F},root::Supernode,p::Array{Int},B::StridedMatrix{F})
    # P A P^T = L L^H and A x = b imply
    #
    #     x = inv(P) (L^H \ (L \ (P b)))
    #
    XNodal = RHSNode{F}(B[p,:],root)
    ForwardSolveLDL!(front,root,XNodal)
    BackwardSolveLDL!(front,root,XNodal)
    XPerm = Unpack(XNodal,root)
    X = XPerm[invperm(p),:]
    X
end
