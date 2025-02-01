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

# Load the graph from a file 
filename_nodes = "PageRank/data/fb-pages-public-figure/fb-pages-public-figure.nodes"
filename_edges = "PageRank/data/fb-pages-public-figure/fb-pages-public-figure.edges"

# Read the nodes
nodes = open(filename_nodes) do file
    readlines(file)
end
# remove the first Line
nodes = nodes[2:end]

# map the nodes to their ids
node_map = Dict()
node_list = []
for (i, node) in enumerate(nodes)
    node = split(node, ",")
    push!(node_list, node[])
end


# Read the edges
edges = open(filename_edges) do file
    readlines(file)
end

# Create a adjacency matrix from the edges 
# nodes come in the form id,name,new_id
# edges come in the form id1,id2

# only take top 1000 nodes



# If your data is zero-indexed, adjust the matrix size accordingly:
n = 1000
A = spzeros(n, n)
# only takes edges with 1...1000



for edge in edges
    edge = split(edge, ",")
    id1 = parse(Int, edge[1]) + 1
    id2 = parse(Int, edge[2]) + 1
    if id1 <= 1000 && id2 <= 1000
        A[id1, id2] = 1
        A[id2, id1] = 1
        println("id1: $id1, id2: $id2")
    end
end
# Check that all nodes have edges
for i in 1:n
    if sum(A[i, :]) == 0
        println("Node $i has no edges")
        # add random edges
        for j in 1:n
            if i != j && rand() < 0.1
                A[i, j] = 1
            end
        end
    end
end

# count edges
edges = 0
for i in 1:n
    edges += sum(A[i, :])
end
println("Number of edges: $edges")



# PageRank iteration

# Initialize the PageRank vector
v = ones(n)/n # teleportation matrix
d = vec(sum(A,dims = 2)) # sum of the columns
dinv = 1 ./ d # inverse of the sum of the columns
Dinv = spdiagm(dinv) # diagonal matrix of the inverse of the sum of the columns
P = A'*Dinv # transition matrix
x = ones(n) / n # PageRank vector
rdata = []
labels = []
for a in [0.1, 0.5, 0.9]
    x = ones(n)/n
    rd = []
    for iteration in 1:20
        rd = [rd ; norm((I - a * P) * x - (1 - a) * v, 2)]
        x = a * P * x + (1-a) * v
    end
    push!(rdata, rd)
    push!(labels, "alpha = $a")
end

plt = plot(title="Residual Norm vs. Iteration", xlabel="Iteration", ylabel="Residual Norm", legend=:topright)
for (i, alp) in enumerate([0.1, 0.5, 0.9])
    plot!(plt, 1:20, rdata[i], label="alp = $(alp)", marker=:o)
end
display(plt)



# Plot 2 - ||x* - xk|| vs. Iteration

data = []

for alp in [0.1, 0.5, 0.9]
    RHS = (1 - alp) * v # right hand side of the equation
    LHS = (I - alp * P) # left hand side of the equation
    x_star = LHS\RHS # solve the equation 
    xd = []
    x = ones(n) / n
    x_init = ones(n) / n
    for k in 1:20
        xd = [xd; norm(x_star - x, 1)]
        x = (alp * P) * x + (1 - alp) * v
    end
    push!(data, xd)
end
plt = plot(title="||x* - xk|| vs. Interation", xlabel="Iteration", ylabel="Residual Norm", legend=:topright)
for (i, alp) in enumerate([0.1, 0.5, 0.9])
    plot!(plt, 1:20, data[i], label="alpha = $alp", marker=:o)
end
display(plt)
