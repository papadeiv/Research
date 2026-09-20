# Import the necessary packages and local modules
include("inc.jl")

# Model parameters
S0 = 0.035

α = 0.12
β = 790.0
Υ = 100*3.15*1e7

Sn_eq = 0.034912
St_eq = 0.035435
Ss_eq = 0.034427
Sip_eq = 0.034668
Sb_eq = 0.034538

Vn = 3.683*1e16
Vt = 5.418*1e16
Vs = 6.097*1e16
Vip = 1.486*1e17
Vb = 9.925*1e17

Ts = 7.919
T0 = 3.87

C = Vn*Sn_eq + Vt*St_eq + Vip*Sip_eq + Vs*Ss_eq + Vb*Sb_eq

λ = 1.62*1e7
γ = 0.36                                           # Control parameter (location of HB_sub)
μ = 22*1e-8

Kn = 1.762*1e6
Ks = 1.872*1e6

Fn(H) = (0.486 + 0.131*H)*1e6
Ft(H) = (-0.997 + 0.696*H)*1e6

Sn0 = (Sn_eq - S0)*100.0
St0 = (St_eq - S0)*100.0
Ss = (Ss_eq - S0)*100.0 
Sb = (Sb_eq - S0)*100.0 
Sip(Sn, St) = 100.0*(C - (Vn*Sn + Vs*Ss + Vb*Sb + Vt*St)/100.0 - S0*(Vb + Vn + Vt + Vip + Vs))/Vip

q(Sn) = λ*(α*(Ts - T0) + (β/100.0)*(Sn - Ss))/(1.0 + λ*α*μ)

# System parameters
H0 = -0.379                                        # Initial value of the bifurcation parameter
Hf = 0.500                                         # Final value of the bifurcation parameter
ε = 1e-4                                           # Slow timescale
σ = 0.001::Float64                                 # Noise level (additive)
D = (σ^2)/2.0                                      # Diffusion level (additive) 
δt = 5e-2                                          # Timestep

# Dynamical system: drift 
f1_pos(Sn, St, H) = (Υ/Vn)*(q(Sn)*(St/100.0 - Sn/100.0) + Kn*(St/100.0 - Sn/100.0) - Fn(H)*S0)
f1_neg(Sn, St, H) = (Υ/Vn)*(abs(q(Sn))*(Sb/100.0 - Sn/100.0) + Kn*(St/100.0 - Sn/100.0) - Fn(H)*S0)
f1(Sn, St, H) = q(Sn) >= 0 ? f1_pos(Sn, St, H) : f1_neg(Sn, St, H)
f2_pos(Sn, St, H) = (Υ/Vt)*(q(Sn)*(γ*Ss/100.0 + (1 - γ)*Sip(Sn, St)/100.0 - St/100.0) + Ks*(Ss/100.0 - St/100.0) + Kn*(Sn/100.0 - St/100.0) - Ft(H)*S0)
f2_neg(Sn, St, H) = (Υ/Vt)*(abs(q(Sn))*(Sn/100.0 - St/100.0) + Ks*(Ss/100.0 - St/100.0) + Kn*(Sn/100.0 - St/100.0) - Ft(H)*S0)
f2(Sn, St, H) = q(Sn) >= 0 ? f2_pos(Sn, St, H) : f2_neg(Sn, St, H)
f = [f1, f2]

# Dynamical system: shift
Λ(t) = ε

# Dynamical system: diffusion
g(x, y) = σ
η = [g, g]

# Define the main algorithm
function main()
        # Solve the ensemble slow-fast SDEs
        equilibria = get_equilibria(f1, f2, H0, guesses=[[Sn0,St0], [-0.5,0.5], [-0.1,0.6]])
        z0 = [equilibria.stable[1]; H0]
        ensemble = evolve(f, η, Λ, z0, stepsize=δt, endparameter=Hf)

        # Extract the x and y components
        x = (ensemble.state[1])[:, 1]
        y = (ensemble.state[1])[:, 2]

        # Extract the parameter's shift and timesteps
        H = ensemble.parameter
        t = ensemble.time

        # Plot the timeseries and phase space trajectory
        fig = Figure()
        ax1 = Axis(fig[1:3,1:2]) 
        lines!(ax1, x, y, color = (:black,0.35), linewidth = 0.5)
        ax2 = Axis(fig[4,1])
        lines!(ax2, t, x, color = :red, linewidth = 1.0)
        ax3 = Axis(fig[4,2])
        lines!(ax3, t, y, color = :blue, linewidth = 1.0)
        save("../../res/fig/wood/trajectory.png", fig)

        # Plot the bifurcation diagram
        xon = Vector{Vector{Float64}}()
        yon = Vector{Vector{Float64}}()
        xoff = Vector{Vector{Float64}}()
        yoff = Vector{Vector{Float64}}()
        xunst = Vector{Vector{Float64}}()
        yunst = Vector{Vector{Float64}}()
        xsddl = Vector{Vector{Float64}}()
        ysddl = Vector{Vector{Float64}}()
        for n in 1:length(H)
                if n == 1
                        equilibria = get_equilibria(f1, f2, H[n], guesses=[[Sn0,St0], [-0.15,0.6], [-0.2,0.5]])
                        push!(xon, [(equilibria.stable[1])[1], H[n]])
                        push!(yon, [(equilibria.stable[1])[2], H[n]])
                        push!(xoff, [(equilibria.stable[2])[1], H[n]])
                        push!(yoff, [(equilibria.stable[2])[2], H[n]])
                        push!(xunst, [(equilibria.unstable[1])[1], H[n]])
                        push!(yunst, [(equilibria.unstable[1])[2], H[n]])
                else
                        ON_previous = [(xon[end])[1], (yon[end])[1]] 
                        OFF_previous = [(xoff[end])[1], (yoff[end])[1]] 
                        UNST_previous = [(xunst[end])[1], (yunst[end])[1]] 
                        equilibria = get_equilibria(f1, f2, H[n], guesses=[ON_previous, OFF_previous, UNST_previous])
                        if length(equilibria.unstable) < 2 && length(equilibria.stable) > 1
                                println("H = $(H[n]): N. stable = $(length(equilibria.stable)), N.unstable = $(length(equilibria.unstable))")
                                equilibria = get_equilibria(f1, f2, H[n], guesses=[ON_previous, OFF_previous, UNST_previous])
                                push!(xon, [(equilibria.stable[1])[1], H[n]])
                                push!(yon, [(equilibria.stable[1])[2], H[n]])
                                push!(xoff, [(equilibria.stable[2])[1], H[n]])
                                push!(yoff, [(equilibria.stable[2])[2], H[n]])
                                push!(xunst, [(equilibria.unstable[1])[1], H[n]])
                                push!(yunst, [(equilibria.unstable[1])[2], H[n]])
                        elseif length(equilibria.unstable) == 2
                                println("H = $(H[n]): N. stable = $(length(equilibria.stable)), N.unstable = $(length(equilibria.unstable))")
                                equilibria = get_equilibria(f1, f2, H[n], guesses=[ON_previous, OFF_previous, UNST_previous])
                                push!(xoff, [(equilibria.stable[1])[1], H[n]])
                                push!(yoff, [(equilibria.stable[1])[2], H[n]])
                                push!(xsddl, [(equilibria.unstable[1])[1], H[n]])
                                push!(ysddl, [(equilibria.unstable[1])[2], H[n]])
                                push!(xunst, [(equilibria.unstable[2])[1], H[n]])
                                push!(yunst, [(equilibria.unstable[2])[2], H[n]])
                        else
                                equilibria = get_equilibria(f1, f2, H[n], guesses=[ON_previous, OFF_previous, UNST_previous])
                                push!(xoff, [(equilibria.stable[1])[1], H[n]])
                                push!(yoff, [(equilibria.stable[1])[2], H[n]])
                        end
                end

        end
        fig = Figure()
        ax = Axis(fig[1,1])
        lines!(ax, last.(xon), first.(xon), color = :blue, linewidth = 2.0)
        lines!(ax, last.(xoff), first.(xoff), color = :blue, linewidth = 2.0)
        lines!(ax, getindex.(xunst, 2), getindex.(xunst, 1), color = :red, linewidth = 2.0)
        lines!(ax, getindex.(xsddl, 2), getindex.(xsddl, 1), color = :goldenrod, linewidth = 2.0)
        lines!(ax, H, x, color = :black, linewidth = 1.0)
        save("../../res/fig/wood/bifurcation_diagram_x.png", fig)

        fig = Figure()
        ax = Axis(fig[1,1])
        lines!(ax, last.(yon), first.(yon), color = :blue, linewidth = 2.0)
        lines!(ax, last.(yoff), first.(yoff), color = :blue, linewidth = 2.0)
        lines!(ax, getindex.(yunst, 2), getindex.(yunst, 1), color = :red, linewidth = 2.0)
        lines!(ax, getindex.(ysddl, 2), getindex.(ysddl, 1), color = :goldenrod, linewidth = 2.0)
        lines!(ax, H, y, color = :black, linewidth = 1.0)
        save("../../res/fig/wood/bifurcation_diagram_y.png", fig)

        #=
        # Loop over the initial parameter values
        for h0 ∈ H0
                # Plot the nullclines
                Sn_domain = LinRange(-1.0,1.0,1000)
                St_domain = LinRange(-1.0,1.0,1000)
                SN = [f1(Sn,St,h0) for Sn in Sn_domain, St in St_domain]
                ST = [f2(Sn,St,h0) for Sn in Sn_domain, St in St_domain]
                fig = Figure()
                ax = Axis(fig[1,1])
                contour!(ax, Sn_domain, St_domain, SN, levels=[0], color=:red, linewidth=3)
                contour!(ax, Sn_domain, St_domain, ST, levels=[0], color=:blue, linewidth=3)
                save("../../res/fig/wood/nullclines/$glb_idx.png", fig)
                global glb_idx = glb_idx + 1
        end
        =#
end

# Execute the main
main()
