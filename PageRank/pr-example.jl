using LinearAlgebra
using MAT
using SparseArrays
using Plots
include("graphs-nopackages.jl")


filename = "PageRank/data/Karate.mat"
M = matread(filename)
A = M["A"]
A = sparse(A)
xy = M["xy"]




f = display_graph(A,xy)
n = size(A,1)

# A sparse matrix is stored using 3 vectors
A.colptr
A.nzval
A.rowval

# Compute PageRank
alp = 0.1

v = ones(n)/n
d = vec(sum(A,dims = 2))
dinv = 1 ./ d
Dinv = spdiagm(dinv)
P = A'*Dinv

G = (I - alp*P)
y = (1-alp)*v
x = round.(G\y,digits=2)

for i = 1:n
    annotate!(f,xy[i,1],xy[i,2]-1.5, text("$(x[i])", :red, 10))
end
f



## Try eigenvector plot
L = nlap(A)
E = eigen(Matrix(L))
V = E.vectors[:,2:5]
xy = [V[:,2] V[:,3]]
display_graph(A,xy)

## Try a cycle
n = 10
A = cycle_graph(n)
L = nlap(A)
E = eigen(Matrix(L))
V = E.vectors[:,2:3]
xy = [V[:,1] V[:,2]]
display_graph(A,xy)


## Try a clique
n = 6
A = clique(n)
L = lap(A)
E = eigen(Matrix(L))
V = E.vectors[:,2:3]
xy = [V[:,1] V[:,2]]
display_graph(A,xy)


