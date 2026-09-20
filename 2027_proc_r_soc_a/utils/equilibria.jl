"""
Compute the equilibria of a one- and two-dimensional vector field with different signatures.

Author: Davide Papapicco
Affil: U. of Auckland
Date: 21-08-2025
"""

function get_equilibria(f::Function, μ::Float64; domain=[-Inf,Inf])
        # Reduce the parameter dependent dynamics to a 1D scalar function
        F(x) = f(x, μ)

        # Find the zeros of the function in the specified interval
        equilibria = find_zeros(F, domain[1], domain[2])

        # Define empty arrays to separate stable from unstable equilibria
        stable = Float64[]
        unstable = Float64[]

        # Loop over the equilibria 
        for eq in equilibria
                # Compute the Jacobian at the equilibrium
                stability = ForwardDiff.derivative(F,eq)
                # Determine the stability by checking the sign 
                if stability < 0
                        push!(stable, eq)
                else
                        push!(unstable, eq)
                end
        end

        # Return all the equilibria 
        return (
                stable = stable,
                unstable = unstable
               ) 
end

function get_equilibria(f::Function; domain=[-Inf,Inf])
        # Find the zeros of the function in the specified interval
        equilibria = find_zeros(f, domain[1], domain[2])

        # Define empty arrays to separate stable from unstable equilibria
        stable = Float64[]
        unstable = Float64[]

        # Loop over the equilibria 
        for eq in equilibria
                # Compute the Jacobian at the equilibrium
                stability = ForwardDiff.derivative(f,eq)
                # Determine the stability by checking the sign 
                if stability < 0
                        push!(stable, eq)
                else
                        push!(unstable, eq)
                end
        end

        # Return all the equilibria 
        return (
                stable = stable,
                unstable = unstable
               ) 
end

function get_equilibria(f1::Function, f2::Function, μ::Float64; guesses=[[0.0, 0.0], [1.0, 1.0], [-1.0, -1.0], [1.0, -1.0], [-1.0, 1.0]])
        # Define the vector field
        function vector_field!(F, u, μ)
                x, y = u[1], u[2]
                F[1] = f1(x, y, μ) 
                F[2] = f2(x, y, μ) 
        end

        # Define empty arrays to separate stable from unstable equilibria
        stable = Vector{Vector{Float64}}()
        unstable = Vector{Vector{Float64}}()

        # Loop over the guesses
        for guess in guesses
                # Solve the nonlinear 0-problem
                prob = NonlinearSolve.NonlinearProblem(vector_field!, guess, μ)
                sol = NonlinearSolve.solve(prob)
                X = sol[1] 
                Y = sol[2] 

                # Avoid exporting duplicates
                if any(norm([X, Y] - [sx, sy]) < 1e-6 for (sx, sy) in vcat(stable, unstable))
                        continue
                end

                # Compute the linearisation (Jacobian)
                J = ForwardDiff.jacobian(u -> [f1(u[1], u[2], μ), f2(u[1], u[2], μ)], [X, Y])
                # Compute the spectrum of the linear field (eigenvalues of the Jacobian)
                eigvals_J = eigen(J).values

                # Establish asymptotic stability
                if all(real.(eigvals_J) .< 0)
                        push!(stable, [X, Y])
                else
                        push!(unstable, [X, Y])
                end
        end

        # Return all the equilibria 
        return (
                stable = stable,
                unstable = unstable
               ) 
end
