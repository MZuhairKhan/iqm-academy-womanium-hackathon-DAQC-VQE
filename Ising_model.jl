using LinearAlgebra, SparseArrays

#=
    Checks if an entry of an a Tuple what represents a point in the cristaline laticce is a valid
        point, if some of entry is not in the defined ranges, the function returns false.
=#
function valid_neighbohrs(point::Tuple{Int64, Int64, Int64},ranges::Vector{UnitRange{Int64}})
    pr = @. (point in ranges)

    pr_t = !(false in pr) # if there is a point with no valid coordinates respect with a certain range, then eliminate

    pr_t
end


#=
    Generates a list of valid neighbohrs positons of certain point in the laticce
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
    Auxiliary function that generates a vector of vector from the raws of a matrix
=#
columns(M) = (view(M, :, i) for i in 1:size(M, 2))


#=
    Auxiliary function that prepares the array for 0 a based index system like Python 
=#
function to_Python(x)
    x .- 1
end


#=
  Calulates from a tuple of dimentions a list of valid connectivities between points in the lattice 
=#
function connectivity(dims::Tuple{Int64,Int64,Int64})
    len=length(dims)

    dimention_range = Array{UnitRange{Int64}}(UndefInitializer(), len)
    
    for i in 1:len
        dimention_range[i] = range(1,dims[i])
    end

    points = [(r1,r2,r3) for r1 in dimention_range[1],r2 in dimention_range[2], r3 in dimention_range[3]]
    
    dic_points = Dict()
    points_dic = Dict()

    i = 1
    for points in points
        push!(dic_points,points=>i)
        push!(points_dic,i=>points)
        i += 1
    end

    points_and_neighbors = Dict()

    for point in points
        index = dic_points[point]
        index_neighbohrs = neighbohrs(point,dimention_range)
        a1 = []#dic_points[index_neighbohrs...]
        for ind_nei in index_neighbohrs
            push!(a1,dic_points[ind_nei])
        end
        push!(points_and_neighbors,index =>a1)
    end
    
    connectivities = []
    
    kys = keys(points_and_neighbors) |> collect
    for ks in kys
        nhbs = points_and_neighbors[ks] 
        for pt in nhbs
            push!(connectivities,[ks,pt])
        end
    end

    connectivities = Matrix(reduce(hcat,connectivities)')

    for i in 1:size(connectivities,1)
        a1, a2 = connectivities[i,:]
        if a1 > a2
            connectivities[i,1] = a2
            connectivities[i,2] = a1
        end
    end

    connectivities = sort(sort(connectivities,dims=2),dims=1)'

    connectivities = columns(connectivities)|> collect |> unique

    to_Python.(connectivities)
end

# Hamiltonian creation function for a Ising model
function create_zz_hamiltonian(num_qubits::Int,connectivity::Vector{Vector{Int}},
    h_coeffs::Vector{Float64})
    
    connectivity = @. (a -> a .+ 1)(connectivity)#[a .+ 1 for a in connectivity]

    dim = 2^num_qubits
    num_connections = length(connectivity)
    zz_hamiltonian = zeros(ComplexF64,dim,dim) |> sparse

    σz = ComplexF64.([1. 0.;0. -1.]) |> sparse

    for c in 1:num_connections
        ops_to_tensor = [1.0*ComplexF64.(I(2)) for i in 1:num_qubits]
        ops_to_tensor[connectivity[c][1]] .= σz
        ops_to_tensor[connectivity[c][2]] .= σz
        ops_to_tensor
        zz_hamiltonian .+= h_coeffs[c] * kron(ops_to_tensor...)
    end

    return zz_hamiltonian
end

#= With no other orguments provided, create_zz_hamiltonian generates a matrix and a number of
qubits (as many as atoms in the lattice), also, supose a homogeneous interaction. We are
using the multiple dispatchs capabilities of Julia. This function receives a tuple with 3 elements and 
returns a hamiltonian Matrix and the number of qbits needed to perform subsequent operations.
=# 
function create_zz_hamiltonian(dims::Tuple{Int64,Int64,Int64})
    
    connectivities= connectivity(dims)
    
    num_qbits = prod(dims)
    h_coeffs = ones(length(connectivities))
    
    zz_hamiltonian = create_zz_hamiltonian(num_qbits,connectivities,h_coeffs)
    
    zz_hamiltonian = real.(zz_hamiltonian)  # We took only the real part
    
    return zz_hamiltonian, num_qbits, h_coeffs, connectivities
end
