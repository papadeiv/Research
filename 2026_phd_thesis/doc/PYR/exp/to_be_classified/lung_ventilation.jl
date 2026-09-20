using DifferentialEquations, LinearAlgebra, LaTeXStrings, ProgressMeter
using PyCall
include("Lattice.jl")

printstyled("LUNG VENTILATION LATTICE MODEL\n"; bold=true, underline=true, color=:light_blue)

# Define the lattice
n = 300
m = 300
N = n*m
L = Lattice(n, m)

# Initial condition of the lattice (low-noise)
#x0 = 0.93::Float64.*ones(Float64, length(L.grid))
#threshold = 1.08
# Initial condition of the lattice (high-noise)
x0 = 0.93::Float64.*ones(Float64, length(L.grid))
threshold = 0.68
# Initial condition of the bifurcation parameter (degree of smooth muscle activation)
K = 0.5::Float64 # smooth muscle activation
push!(x0, K) 
# SDE's parameters
σ = 0.080::Float64 # noise level
push!(x0, σ)
ε = 0.010::Float64 # parameter's drift
push!(x0, ε)

# Model's parameters
k = 14.1::Float64 # (airway) smooth muscle mass
A = 0.63::Float64 # inter-airway coupling
Pi = 0.96::Float64 # inflection point of pressure-radius interdependence
Pb = 0.0::Float64 # breathing pressure
Pb0 = 7.25::Float64 # breating pressure's IC

# Sliding window's parameters
T = 100.00
δt = 1e-2
width = 20

# Definition of the dynamics (lung ventilation) function
sigmoid(x::Real) = one(x)/(one(x) + exp(-x))
# Definition of the SDE - deterministic part
function iip_det!(f, x, neigh, t)
        sum = 0.0::Float64
        for j=1:N
                sum += x[j]^4
        end
        Pb = (Pb0*N)/(sum)
        for j=1:N
                r = -Pb + x[N+1]*(k/x[j]) - Pb*A*(x[j]^4 + x[neigh[1,j]]^4 + x[neigh[2,j]]^4 + x[neigh[3,j]]^4 + x[neigh[4,j]]^4)*(1.0::Float64 - x[j] + 1.5::Float64*(1-x[j])^2) + Pi
                f[j] = sigmoid(-r) - x[j]
        end
        f[N+1] = ε
        return nothing
end
# Definition of the SDE - stochastic part
function iip_stoc!(f, x, neigh, t)
        for j=1:N
                f[j] = σ
        end
        f[N+1] = 0.0::Float64
        return nothing
end
problem = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T), L.connectivity)

# Compute the forward-time trajectory
solution = solve(problem, EM(), dt=δt, saveat=0.1)
time = solution.t
states = solution.u

# Compute the spatial aggregate of the solutions
aggregate = [(1/N)*sum(states[t][1:end-3]) for t=1:length(time)] 
# Extract the parameter values at each timestep
μ = [states[t][end-2] for t=1:length(time)]

# Plot the timeseries
CairoMakie.activate!(; px_per_unit = 3)
fig = Figure(size = (1800, 600), backgroundcolor = :transparent)
ax = Axis(fig[1,1],
                # Background
                backgroundcolor = :white,
                spinewidth = 5.0,
                xtickwidth = 5.0,
                ytickwidth = 5.0,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((μ[1],μ[end]), (aggregate[end],1)),
                # Title
                #title = (L"\mu=%$(round(Y; digits=3))"),
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # x-axis
                xlabel = L"\mathbf{\mu}",
                xlabelvisible = true,
                xlabelsize = 60,
                xlabelcolor = :black,
                xlabelpadding = -60.0,
                xticks = [0.5, μ[end]],
                xticksvisible = true,
                xticksize = 20,
                xticklabelsvisible = true,
                xticklabelsize = 50,
                xtickformat = "{:.1f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\textbf{mean field}",
                ylabelvisible = true,
                ylabelsize = 60,
                ylabelcolor = :black,
                ylabelpadding = -60.0,
                yticks = [aggregate[end],1],
                yticksvisible = true,
                yticksize = 20,
                xticklabelalign = (:right, :top),
                yticklabelsvisible = true,
                yticklabelsize = 50,
                ytickformat = "{:.1f}",
                yscale = identity,
                yaxisposition = :left,
)
lines!(ax, μ, aggregate, linewidth = 3.5, color = :black)
# High-noise threshold
#lines!(ax, threshold.*ones(2), [aggregate[argmin(aggregate)],1], color = :black, linestyle = :dash)
save("../results/lattice/lung_ventilation/fig3.1.png", fig)

# Postprocessing of the EWS using Python's timeseries classes 
#=
using PyCall
py"""
import sys
#print(sys.version)
sys.path.append('../timeseries')

from Process import Process
from Estimator import Estimator 
from TimeSeries import TimeSeries
from matplotlib import pyplot as plt

import numpy as np

def analyse(data):
        ts = Process()
        ts.realizations = TimeSeries(realizations=data) 
        np.savetxt('./aggregate.csv', ts.realizations.ts, delimiter=',')

        ts.detrend(mode='EMD', order=10)
        imfs = ts.imfs

        np.savetxt('./imfs1.csv', imfs[:,(imfs.shape)[1]-1], delimiter=',')
        np.savetxt('./imfs2.csv', imfs[:,(imfs.shape)[1]-2], delimiter=',')
        np.savetxt('./imfs3.csv', imfs[:,(imfs.shape)[1]-3], delimiter=',')
        np.savetxt('./imfs4.csv', imfs[:,(imfs.shape)[1]-4], delimiter=',')
        np.savetxt('./imfs5.csv', imfs[:,(imfs.shape)[1]-5], delimiter=',')
        np.savetxt('./imfs6.csv', imfs[:,(imfs.shape)[1]-6], delimiter=',')
        np.savetxt('./imfs7.csv', imfs[:,(imfs.shape)[1]-7], delimiter=',')

        EWS1 = Estimator(ts.detrended.ts)
        EWS1.variance(10)
        np.savetxt('./emd.csv', EWS1.var, delimiter=',')

        EWS2 = Estimator(ts.realizations.ts)
        EWS2.variance(10)
        np.savetxt('./nodetrend.csv', EWS2.var, delimiter=',')

        fig = plt.figure(figsize=[12.8,9.6], dpi=200, layout='tight')
        if imfs.shape[1]<=4:
                for n in range(imfs.shape[1]):
                        ax = plt.subplot2grid((4,4), (n,0), colspan=4, fig=fig)
                        ax.plot(np.linspace(0,100,imfs.shape[0]), imfs[:,0], label='IMF 1')

        elif imfs.shape[1]<=6:
                for n in range(imfs.shape[1]):
                        quotient, remainder = divmod(n, 3)
                        ax = plt.subplot2grid((3,2), (remainder,quotient), colspan=1, fig=fig)
                        ax.plot(np.linspace(0,100,imfs.shape[0]), imfs[:,n], label='IMF 1')

        else:
                for n in range(imfs.shape[1]):
                        quotient, remainder = divmod(n, 4)
                        ax = plt.subplot2grid((4,4), (remainder,2*quotient), colspan=2, fig=fig)
                        ax.plot(np.linspace(0,100,imfs.shape[0]), imfs[:,n], label='IMF 1')


        plt.savefig("../results/lattice/lung_ventilation/IMFs.png", transparent=True, dpi=600)
"""
# Analyse and detrend timeseries data
pyanalyse = py"analyse"
pyanalyse(aggregate)
=#

# Dynamic Mode Decomposition
EWS = zeros(Float64, (length(time)-width))
Mode = Array{Float64,2}(undef, N, (length(time)-width))
#=
println("Performing the analysis of the DMD modes")
@showprogress for j=1:(length(time)-width)
        # Assemble the first snapshot matrix
        idx = 1::Int
        X = Array{Float64,2}(undef, N, width)
        for k=j:(j-1)+width
                X[:,idx] = states[k][1:N]
                idx += 1
        end
        # Perform the DMD
        F = svd(X)
        # Assemble the second snapshot matrix
        idx = 1::Int
        Y = Array{Float64,2}(undef, N, width)
        for k=(j+1):j+width
                Y[:,idx] = states[k][1:N]
                idx += 1
        end
        # Perform the eigendecomposition
        S = (F.U)'*Y*(F.Vt)'*inv(diagm(F.S)) 
        Λ = eigvals(S)
        Q = eigvecs(S)
        # Get leading eigenvalue and eigenvector as spatial EWS
        EWS[j] = (findmax(real(Λ)))[1]
        Mode[:,j] = F.U*real(Q[:,(findmax(real(Λ)))[2]])
end
=#

# Export images of the lattice time evolution
show(L, states, EWS, Mode)
