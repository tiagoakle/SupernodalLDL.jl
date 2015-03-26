# SupernodalLDL

Implementation of a supernodal factorization routine for sparse quasi-definite matrices. 
The Cholesky version of this code was developed by Jack Poulson for his 
Advanced Topics in Numerical Linear Algebra Course at ICME Stanford university, Winter 2015.
http://web.stanford.edu/~poulson/courses/ATiNLA15/

We have gathered the code, and formed a self-contained package. We have added comments and 
made minor modifications. 

Using Nick Henderson's QuasiDefinite.jl bindings for his Fortran extension to LAPACK for dense
quasi-definite matrix factorization, we strive to create a Supernodal LDL factorization for 
sparse quasi-definite matrices.

## Clone and Build in Julia

```julia
julia> Pkg.clone("https://github.com/tiagoakle/SupernodalLDL.jl.git")
```

## Load in Julia

```julia
julia> using SupernodalLDL
```
