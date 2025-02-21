#=
PageRank Algorithm:
- Random Surfer Model:
    - \alpha = (0, 1) - following uniform random link
    - 1 - \alpha - teleportation
=#
using LinearAlgebra
using MAT
using SparseArrays
using Plots
include("graphs-nopackages.jl")


# Randomly generate a matrix representation of a graph
function generate_graph(n)
    A = zeros(n,n)
    for i = 1:n
        for j = 1:n
            if i != j
                A[i,j] = rand(0:5)
                if A[i,j] == 1
                    A[i,j] = 1
                else
                    A[i,j] = 0
                end
            end
        end
    end
    return A
end

# Auto gen graph
A = generate_graph(20)
A = sparse(A)
# generate x y coordinates for the nodes from random numbers from 1..100
xy = rand(20,2) * 100
n = 20


A.colptr
A.nzval
A.rowval

# PageRank
alp = 0.9 # alpha value 
v = ones(n)/n # teleportation matrix
d = vec(sum(A,dims = 2)) # sum of the columns
dinv = 1 ./ d # inverse of the sum of the columns
Dinv = spdiagm(dinv) # diagonal matrix of the inverse of the sum of the columns
P = A'*Dinv # transition matrix

v_et =  v * ones(1,n) # teleportation matrix
M = alp * P + (1 - alp) * v_et # matrix M

RHS = (1 - alp) * v # right hand side of the equation
LHS = (I - alp * P) # left hand side of the equation
x = round.(LHS\RHS,digits=2) # solve the equation 
f = display_graph(A,xy)

for i = 1:n
    annotate!(f,xy[i,1],xy[i,2]-1.5, text("$(x[i])", :red, 10))
end
f