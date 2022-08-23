#=
This code generalizes the generation of Issing Hamiltonians given a set of
connectivities and/or a set of connectivities coefficients. It can generate
Hamiltonians up to 12 atoms, or 12 qubits. For example, if we have a line of atoms

A-A-A-A

with positions 

0-1-2-3

their connectivity set, without repeated interactions should be

[[0,1],[1,2],[3,4]]

Another example, if we have a certain connectivity in 1D, for a star-like set as

          A(5)
          |
A(1)-A(2)-A(0)-A(3)-A(4) 
          |
          A(6)

The set of connectivities to be provided to the functions are

[[1,2],[3,4],[0,2],[0,3],[0,6],[0,5]]


Also if we consider a simple cubic lattice (P), and a cristal of 3D dimensions (x,y,z) up to 12 atoms, 
the code can resolve how many qubits do we need to simulate the hamiltonian, which would be 
the unique connectivities, and a set of uniform interaction coefficients.
=#

using LinearAlgebra, SparseArrays # Packages, both Julia native. 

# In Julia isn't mandatory to set the type of the entries of the function, but it improves 
# the speed of compilation and execution of the functions. Also helps to debugs 
# only certain compiled methods and not the whole function, taking advantage of the 
# multiple dispatch capabilities of Julia.

#=
    Checks if an entry of an a Tuple what represents a point in the cristaline laticce is a valid
        point, if some of entry is not in the defined ranges, the function returns false.
    Receives a 3D tuple as a point an a set of ranges. 
=#
function valid_neighbohrs(point::Tuple{Int64, Int64, Int64},ranges::Vector{UnitRange{Int64}}) 
    # Its not mandatory in Julia to set the variable type in the entries of a functions, but
    # doing that we ensure a better process of debugging and shorter compilation and execution times. 
    pr = @. (point in ranges)  # Checks if every entry resides between the ranges

    pr_t = !(false in pr) # if there is a point with no valid coordinates respect with a certain range,
                          # it will be signaled with a false signal

    pr_t      # In Julia we are not force to use the clause "return" if we want to return the last
              # line of the function's code. 
end


#=
    Generates a list of valid neighbohrs positons of certain point in the laticce.
    Receives a 3D tuple point, and a set of ranges. Returns an array of valid 
    neighbohrs (tuples) given a point.
=#
function neighbohrs(point::Tuple{Int64,Int64,Int64},ranges::Vector{UnitRange{Int64}})
    v1 = @. point + (1,0,0)
    v2 = @. point - (1,0,0)
    v3 = @. point + (0,1,0)
    v4 = @. point - (0,1,0)
    v5 = @. point + (0,0,1)
    v6 = @. point - (0,0,1)

    va1 = [v1,v2,v3,v4,v5,v6]
    va1t= [valid_neighbohrs(v,ranges) for v in va1]
    va1[va1t]
end


#=
    Auxiliary function that generates a vector of vector from the raws of a matrix.
        Takes colum slices from a matrix
=#
columns(M) = (view(M, :, i) for i in 1:size(M, 2))


#=
    Auxiliary function that prepares the array for 0 a based index system like Python 
=#
function to_Python(x)
    x .- 1
end


#=
  Calulates from a tuple of dimensions a list of valid connectivities between points in the lattice 
=#
function connectivity(dims::Tuple{Int64,Int64,Int64})
    len=length(dims)

    # Array initialization
    dimension_range = Array{UnitRange{Int64}}(UndefInitializer(), len)
    
    # Generate a vector of ranges, e.g. like [1:2,1:2,1:3] for a crystal lattice with 12 atoms
    for i in 1:len
        dimension_range[i] = range(1,dims[i])
    end

    # generate the list of all possible positons in the lattice
    points = [(r1,r2,r3) for r1 in dimension_range[1],r2 in dimension_range[2], r3 in dimension_range[3]]
    
    # Initialization of dictionaries to set integer names to the atoms positions, eg. 1 => (1,1,1)
    dic_points = Dict()
    points_dic = Dict()

    i = 1
    for points in points # We number the points inside the dictionaries
        push!(dic_points,points=>i)  #Points to number (x,y,z) => #
        push!(points_dic,i=>points)  #Number to point # => (x,y,z)
        i += 1
    end

    points_and_neighbors = Dict()  # A dictionary with a point and his neighbohrs

    for point in points
        index = dic_points[point]  # We get the index of a point 
        index_neighbohrs = neighbohrs(point,dimension_range) # We get the array of valid points
                                                             # wigen the array of ranges and a central point
        a1 = []  # We populate an array with the integer indices of the neighbohrs
        for ind_nei in index_neighbohrs
            push!(a1,dic_points[ind_nei])
        end
        push!(points_and_neighbors,index =>a1)
    end
    
    connectivities = []
    
    kys = keys(points_and_neighbors) |> collect  # We resolve all the different integer labels for our atoms
    for ks in kys # for each labeled atom, we create an list of all connectivities, even the repeated
        nhbs = points_and_neighbors[ks] 
        for pt in nhbs
            push!(connectivities,[ks,pt])
        end
    end

    connectivities = Matrix(reduce(hcat,connectivities)')  # Reduction of the matrix to a 1D vector of vectors

    for i in 1:size(connectivities,1)   # Arrangement of the entries, so we can compare the connectivities.
        a1, a2 = connectivities[i,:]    # so if there are [1,2] and [2,1] in the array, we rearrange the last
        if a1 > a2                      # connectivity to [1,2], so the array nows have 2 equal elements.
            connectivities[i,1] = a2
            connectivities[i,2] = a1
        end
    end

    connectivities = sort(sort(connectivities,dims=2),dims=1)'  # We sort the connectivities araray

    connectivities = columns(connectivities)|> collect |> unique  # We filter out the repeated connectivities

    to_Python.(connectivities)  # Changes the indices to be compatible with Python. 
end

# Hamiltonian creation function for a Ising model
# It receives the number of qubits, a vector of vectors of a set of valid connectivities,
# and a vector of connectivity coefficients.
function create_zz_hamiltonian(num_qubits::Int,connectivity::Vector{Vector{Int}},
    h_coeffs::Vector{Float64})
    
    connectivity = @. (a -> a .+ 1)(connectivity)#[a .+ 1 for a in connectivity] #This code is equivalent

    dim = 2^num_qubits   # Set the number of dimensions
    num_connections = length(connectivity) # Measures the needed number of connections 
    zz_hamiltonian = zeros(ComplexF64,dim,dim) |> sparse # Initializes the Hamiltonian matrix, as a complex sparse matrix
                                                         # A sparse matrix is faster to compute as almost everything is zero.

    σz = ComplexF64.([1. 0.;0. -1.]) |> sparse

    for c in 1:num_connections  # Calculates the terms of the hamiltonian
        ops_to_tensor = [1.0*ComplexF64.(I(2)) for i in 1:num_qubits]
        ops_to_tensor[connectivity[c][1]] .= σz
        ops_to_tensor[connectivity[c][2]] .= σz
        ops_to_tensor
        zz_hamiltonian .+= h_coeffs[c] * kron(ops_to_tensor...)  # we concatenate all the matrices in the ops_to_tensor array and make a kroneker 
                                                                 # vectorial product.
    end

    return zz_hamiltonian
end

#= With no other orguments provided, create_zz_hamiltonian generates a matrix and a number of
qubits (as many as atoms in the lattice), also, supose a homogeneous interaction. We are
using the multiple dispatchs capabilities of Julia. This function receives a tuple with 3 elements and 
returns a hamiltonian Matrix, number of qubits, interaction coefficients, and a vector of connectivities. 
=# 
function create_zz_hamiltonian(dims::Tuple{Int64,Int64,Int64})
    
    connectivities= connectivity(dims)
    
    num_qbits = prod(dims)
    h_coeffs = ones(length(connectivities))
    
    zz_hamiltonian = create_zz_hamiltonian(num_qbits,connectivities,h_coeffs)
    
    return zz_hamiltonian, num_qbits, h_coeffs, connectivities
end


#= Given a vector of connectivities (with Python indexing [starting with 0]), generate a hogeneous Ising hamiltonian automatically.
returns a hamiltonian Matrix, number of qubits, interaction coefficients, and the vector of connectivities. =#
function create_zz_hamiltonian(connectivities::Vector{Vector{Int}})
    num_qubits = length(unique(vcat(connectivities...)))
    h_coeffs = ones(num_qubits)
    
    zz_hamiltonian = create_zz_hamiltonian(num_qubits::Int,connectivities::Vector{Vector{Int}},
    h_coeffs::Vector{Float64})
    return zz_hamiltonian, num_qubits, h_coeffs, connectivities
end