# Import packages
using LinearAlgebra, BandedMatrices

# Include functions
include("../src/build_neighbours.jl")

# Define structs
struct Lattice 
        rows::Int
        cols::Int
        grid::Array{Int,2}
        connectivity::Array{Int,2}
        Lattice(rows, cols) = new(rows, cols, LinearIndices((1:rows,1:cols)), build_neighbours((rows,cols)))
end

include("../src/get_neighbours.jl")
