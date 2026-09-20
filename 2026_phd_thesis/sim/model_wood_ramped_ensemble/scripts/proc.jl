"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Arbitrary cubic potential
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)

# Data structures
solutions = Vector{Vector{Float64}}()           # Solutions of the inference method 
results = Vector{Float64}()                     # Statistical distribution over the ensemble
glb_idx = 1::Integer                            # Global index to loop over the parameter values

function analyse(solutions)
        # Compute estimated stable and unstable equilibria of the cubic
        xs = +(1/(3*solutions[3]))*(sqrt((solutions[2])^2 - 3*solutions[1]*solutions[3]) - solutions[2])
        xu = -(1/(3*solutions[3]))*(sqrt((solutions[2])^2 - 3*solutions[1]*solutions[3]) + solutions[2])

        # Compute the ews
        ΔV = abs(V(xu, solutions) - V(xs, solutions))
        push!(results, exp(-ΔV))
end

function export_data()
        # Extract the coefficients 
        mat = transpose(reduce(hcat, solutions))
        μ, c1, c2, c3 = eachcol(mat)

        # Construct the dataframe with its header
        df = DataFrame(parameter = μ, c1 = c1, c2 = c2, c3 = c3, LDP = results)

        # Export the matrix in csv
        CSV.write("../../res/data/wood/$glb_idx.csv", df)
end

function import_data(idx)
        # Import the least-squares solution of the potential reconstruction
        df = CSV.read("../../res/data/wood/$idx.csv", DataFrame)
        μc = df.parameter
        c1 = df.c1
        c2 = df.c2
        c3 = df.c3
        LDP = df.LDP
      
        return (
                coefficients = [[μc[n], c1[n], c2[n], c3[n]] for n in eachindex(μc)],
                ews = LDP 
               )
end
