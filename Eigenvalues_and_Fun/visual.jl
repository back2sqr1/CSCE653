using LinearAlgebra
using MAT
using SparseArrays
using Plots
using MultivariateStats
include("graphs-nopackages.jl")
using LinearAlgebra

filename = "Eigenvalues_and_Fun/dolphins.mat"
M = matread(filename)

# print the keys of the dictionary
println(keys(M["Problem"]))
A = M["Problem"]["A"]

# get xy values for the nodes
g = Graph(A)
gplot(g, layout=spring_layout)

norm_lap = Matrix(nlap(A))


# get second and third Eigenvalues 
eigvecs(norm_lap)
eigenvalues,eigenvectors = eigen(norm_lap)

fiedler_value = eigenvalues[60]
fiedler_vector = eigenvectors[:,60]

# fiedler cut red(positive) and blue(negative)
nodefillc = []
for i in 1:nv(g)
    if fiedler_vector[i] > 0
        push!(nodefillc, colorant"red")
    else
        push!(nodefillc, colorant"blue")
    end
end

gplot(g, nodefillc=nodefillc, layout=spring_layout)