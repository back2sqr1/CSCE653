using LinearAlgebra
using MAT
using SparseArrays
using Plots
include("graphs-nopackages.jl")
using LinearAlgebra
using Graphs
using GraphPlot


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

fiedler = eigenvalues[2]
fiedler_vector = eigenvectors[:,2]

# degree matrix
d = vec(sum(A,dims = 2))
D = Diagonal(d)

indices = collect(1:length(fiedler_vector))
y = fiedler_vector ./ sqrt.(d)

# the y values and keep track of indicies
sorted_indices = sortperm(y)
y = sort(y)



# sweep cut through nodes
set = []
mn_conductance = Inf
best_set = []
for i in sorted_indices
    set = append!(set, i)
    complement = setdiff(indices, set)
    # print("Set: $set Complement: $complement\n")
    cut_cnt = 0
    for j in set
        for k in complement
            if A[j,k] == 1
                cut_cnt += 1
            end
        end
    end
    vol_set = 0
    for j in set
        for k in set 
            if j != k
                vol_set += A[j,k]
            end
        end
    end
    vol_complement = 0
    for j in complement
        for k in complement
            if j != k
                vol_complement += A[j,k]
            end
        end
    end
    if vol_complement == 0 || vol_set == 0
        continue
    end
    conductance = cut_cnt / min(vol_set, vol_complement)
    # println("Conductance Candidate: $conductance $vol_complement $vol_set")
    if conductance < mn_conductance
        mn_conductance = conductance
        best_set = copy(set)
    end
end
mn_conductance
println("Best set: $best_set")

# color the nodes in the best set red and blue otherwise
nodefillc = []
for i in 1:nv(g)
    if i in best_set
        push!(nodefillc, colorant"red")
    else
        push!(nodefillc, colorant"blue")
    end
end

gplot(g, layout=spring_layout, nodefillc=nodefillc)