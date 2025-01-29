"""
Given an adjacency matrix and a 2-D layout, plot the graph
"""
function display_graph(A::SparseMatrixCSC{Float64,Int64},xy::Matrix{Float64},grayscale = 0.0,ms = 6,lw = 1)
    f = plot(legend=false, axis = false,grid = false,xticks = false,yticks = false)
    ei,ej,w = findnz(triu(A))
    scatter!(f,xy[:,1],xy[:,2],color = RGB(grayscale,grayscale,grayscale),markersize = ms, markerstrokecolor =  RGB(grayscale,grayscale,grayscale))
    lx = [xy[ei,1]';xy[ej,1]']
    ly = [xy[ei,2]';xy[ej,2]']
    for i = 1:length(w)
        # draws line from the first point, to the second point
        plot!(f,lx[:,i],ly[:,i],color = RGB(grayscale,grayscale,grayscale), linewidth = lw)
    end
    return f
  end

  """
Given an adjacency matrix, add an edge (i,j) with weight v.
If there alreay was an edge, this will overwrite it.

This assumes the graph is undirected, so the matrix should
be symmetric.

The exclamation point symbol in julia function is the convention for indicating
    that this changes the array that it passed in.

The "= 1.0" indicates that the default value for v is 1.0.
"""
function change_edge!(A,i,j,v = 1.0)
    A[i,j] = v
    A[j,i] = v
end


"""
Read in a matrix A in in matrix market format:

First line is M, N, L (rows, columns, and number of nonzeros)
Subsequent lines are of the format
    i, j, v (row index, column index, v = value of A(i,j))
"""
function read_mm_matrix(filename)

    f = readlines(filename)  # reads in the file, f[i] will be line if

    ind = 1
    searching_for_start = true

    # iterate past the comments parts
    while searching_for_start
        if f[ind][1] == '%'
            ind += 1
        else
            searching_for_start = false
        end
    end
    @show ind
    mnl = f[ind]              # this gives the first line that isn't a comment
    mnl = split(mnl)          # splits based on spaces so each number is an entry
    mnl = parse.(Int64,mnl)   # parses the strings into integers. The . notation means "do this to each entry"

    # allocate space for storing (i,j,v) information
    I = zeros(Int64,mnl[3])
    J = zeros(Int64,mnl[3])
    V = zeros(Int64,mnl[3])
    for k = ind+1:length(f)
        entry = parse.(Int64,split(f[k]))
        I[k-ind] = entry[1]
        J[k-ind] = entry[2]
        if length(entry) == 3
            V[k-ind] = entry[3]
        else
            V[k-ind] = 1.0
        end
    end


    # check that number of rows equals number of columns.
    # We want this because we are assuming this is an adjacency matrix for a graph
    @assert(mnl[1] == mnl[2])
    n = mnl[1] # number of nodes
    A = sparse(I,J,V,n,n)
    A = A+A'
    return A
end

function cycle_graph(n)
    A = sparse(collect(1:n),[collect(2:n);1],ones(n),n,n)
    return max.(A,A')
end


function clique(n)
    A = ones(n,n)
    for i = 1:n
        A[i,i] = 0
    end
    return sparse(A)
end

"""
From: https://github.com/dgleich/flow-and-eigenvectors
"""
function simple_spectral_eigenvectors(A,k;nsteps=500,tol=1e-6)
  @assert issymmetric(A) # expensive, but useful...
  n = size(A,1)
  d = vec(sum(A,dims=1))
  nd2 = d'*d
  X = randn(n,k)
  # project x orthogonal to d
  X .-= ((d'*X)/nd2).*d
  #x ./= x'*(d.*x) # normalize
  Q = qr!(sqrt.(d).*X).Q
  X = sqrt.(1.0./d).*Matrix(Q)
  for i=1:nsteps
    X .+= (A*X)./d     # matrix-vector with (I + AD^{-1}) X
    X .-= ((d'*X)/nd2).*d
    Q = qr!(sqrt.(d).*X).Q
    X = sqrt.(1.0./d).*Matrix(Q)
  end
  # make sure the first component is positive
  X .*= repeat(sign.(X[1,:])',size(X,1),1)

  return X
end


function nlap(A)
    d = vec(sum(A,dims = 2))
      Dhalf = Diagonal(d.^(-1/2))
      L = I - Dhalf*A*Dhalf
  end
  
  function lap(A)
    d = vec(sum(A,dims = 2))
    D = Diagonal(d)
    L = D - A
  end


function dirgraphplot(A,xy)

    # Extract x and y coordinates
    node_coords = xy
    adj_matrix = A
    x = node_coords[:, 1]
    y = node_coords[:, 2]

    # Create the plot
    f = plot(x, y, xlim = (minimum(x) - 1/2, maximum(x) + 1/2), ylim = (minimum(y) - 1/2, maximum(y) + 1/2), aspect_ratio = 1, seriestype = :scatter, markersize = 15, markercolor = :white,label = "", legend = false,axis = false,grid = false)

    # Add edges with arrows
    for i in 1:size(adj_matrix, 1)
        for j in 1:size(adj_matrix, 2)
            if adj_matrix[i, j] == 1
                # Start and end points
                x_start, y_start = x[i], y[i]
                x_end, y_end = x[j], y[j]
                
                # Draw an arrow
                arrow_x = [x_start, x_end]
                arrow_y = [y_start, y_end]
                plot!(f, arrow_x, arrow_y, arrow = :arrow, lw = 2, label = "",color = :black)
            end
        end
    end
    annotate!(f,[(x[i], y[i], text("$(i)", :center,:red)) for i in 1:length(x)])

    return f
end


function dirgraphplot(A,xy,S)

    # Extract x and y coordinates
    node_coords = xy
    adj_matrix = A
    x = node_coords[:, 1]
    y = node_coords[:, 2]

    # Create the plot
    f = plot(x, y, xlim = (minimum(x) - 1/2, maximum(x) + 1/2), ylim = (minimum(y) - 1/2, maximum(y) + 1/2), aspect_ratio = 1, seriestype = :scatter, markersize = 15, markercolor = :white,label = "", legend = false,axis = false,grid = false)

    # Add edges with arrows
    for i in 1:size(adj_matrix, 1)
        for j in 1:size(adj_matrix, 2)
            if adj_matrix[i, j] == 1
                # Start and end points
                x_start, y_start = x[i], y[i]
                x_end, y_end = x[j], y[j]
                
                # Draw an arrow
                arrow_x = [x_start, x_end]
                arrow_y = [y_start, y_end]
                plot!(f, arrow_x, arrow_y, arrow = :arrow, lw = 2, label = "",color = :black)
            end
        end
    end
    annotate!(f,[(x[i], y[i], text("$(i)", :center,:black)) for i in 1:length(x)])
    plot!(f,x[S], y[S], seriestype = :scatter, markersize = 15, markercolor = :blue,label = "", legend = false)
    return f
end