import Base: SparseMatrixCSC
using Base.Test
include("./SnodeLDL.jl")

function laplacian2d(nx,ny)
    M = speye(nx*ny)
    for x=1:nx
        for y=1:ny
            s = x+(y-1)*nx
            M[s,s] = 4
            if x > 1 
                M[s,s-1] = -1
            end
            if x < nx 
                M[s,s+1] = -1 
            end
            if y > 1
                M[s,s-nx] = -1 
            end
            if y < ny 
                M[s,s+nx] = -1
            end
        end
    end
    M
end
# The grid dimensions for the 2D finite-difference mesh
nx=ny=40;
n=nx*ny;
# The maximum acceptable leaf-node size
leafSize = 64;
A = laplacian2d(nx,ny);
pND, pND_inv = nodeND(A);

n = A.m;
g = sparseToAdj(A);

# Ensure that Supernode is compiled before timing it
nSmall = 10;
ASmall = speye(nSmall,nSmall);
gSmall = sparseToAdj(ASmall);
pSmall = [1:nSmall];
rootSmall = Supernode(gSmall,pSmall);

p = [1:n];
@time root = Supernode(g,p);
APerm = A[p,p];

# Ensure that Front{Float64} is compiled
frontSmall = Front{Float64}(ASmall[pSmall,pSmall],rootSmall);

# Time the construction after compilation has been guaranteed
@time front = Front{Float64}(APerm,root);

# Run Cholesky on the small example to ensure that it is compiled
Cholesky!(frontSmall);

# Time the actual sparse-direct Cholesky factorization
@time Cholesky!(front);

# Run Solve on a small example to ensure that it is compiled
BSmall = randn(nSmall,10);
Solve(frontSmall,rootSmall,pSmall,BSmall);

# Now time the actual instance of the solve
B = randn(n,10);
@time X = Solve(front,root,p,B);
println( "Relative residual norm $(norm(B-A*X)/vecnorm(B))) ");
