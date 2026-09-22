include("../../../../inc/PlottingTools.jl")
include("../../../../inc/PotentialLearning.jl")

# Define the domain of the observations
a = -1.000 
b = -a 
domain = [a, b]

# Define the deterministic function
f(x) = x^2

# Number of samples
N_samples = convert(Int64,1e1)

# Number of standard deviations to be simulated
Nσ = convert(Int64,2e0)

# Vector of different standard deviations of the samples
σ = LinRange(0.100, 0.000, Nσ)

# Define an array to store the residuals
residual = Vector{Float64}(undef, Nσ)

# Loop over the standard deviations
using ProgressMeter
@showprogress for n in 1:Nσ
        # Generate random samples
        z = σ[n].*randn(N_samples)
        
        # Generate noisy function 
        x = LinRange(domain[1], domain[end], N_samples)
        y = [f(xn) for xn in x] + z 

        # Fit a parabola using Polynomial's built-in function
        V = approximate_potential(x, y, degree=2) 
        #display(V.coeffs)
        # Fit a parabola using the implemented least-squares fit
        U, r = fit_potential(x, y, degree=2)
        #display(U.coeffs)
        #display(std(z))
        #display(norm(r))
        #display(mean(abs.(r)))

        # Compute the error function between the least-squares solution and the analytical function
        e(x) = f(x) - U(x)

        # Compute the residual (L2-norm of the error function)
        N_points = 1001
        dx = (domain[end]-domain[1])/(N_points-1)
        residual[n] = sqrt(dx*sum([(e(x))^2 for x in LinRange(domain[1], domain[2], N_points)]))

        # Create the figure for th samples, the analytcial function and the least-square solution
        local fig, ax = mkfig()
        set_ticks(ax, x, y, n_ticks = 3)

        # Plot the Polynomial's fit
        lines!(ax, LinRange(domain[1], domain[end], 1000), [V(x) for x in LinRange(domain[1], domain[end], 1000)], linewidth = 5, color = (:black, 0.3))
        # Plot the implemented least-squares fit
        lines!(ax, LinRange(domain[1], domain[end], 1000), [U(x) for x in LinRange(domain[1], domain[end], 1000)], linewidth = 3, color = :blue, linestyle = :dash)
        # Plot the analytical function
        lines!(ax, LinRange(domain[1], domain[end], 1000), [f(x) for x in LinRange(domain[1], domain[end], 1000)], linewidth = 3, color = :red)
        # Plot the error function between the analytical one and the least-squares solution
        lines!(ax, LinRange(domain[1], domain[end], 1000), [e(x) for x in LinRange(domain[1], domain[end], 1000)], linewidth = 3, color = :black)
        # Plot the samples
        scatter!(ax, x, y, markersize = 20, color = :blue, strokewidth = 3, strokecolor = :black)
        # Plot the vertical lines from the samples to the least-square solution using the residual
        for m in 1:N_samples
                lines!(ax, [x[m], x[m]], [y[m], y[m]-r[m]], linewidth = 3, color = :blue, linestyle = :dash)
        end

        # Export the figure
        save("../fig/$n.png", fig)
end

# Create the figure for the decay of the residual
fig, ax = mkfig(lab = [L"\sigma", L"||r(x)||_2"],
                ax_orientation = [true, false],
                ticks_lab_trunc = [2,1])
set_ticks(ax, σ, residual, n_ticks = 3)

# Plot the residual decay
lines!(ax, σ, residual, linewidth = 5, color = :red)

# Export the figure
save("../fig/decay.png", fig)
