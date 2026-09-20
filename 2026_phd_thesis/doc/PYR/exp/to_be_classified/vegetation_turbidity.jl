using DifferentialEquations, LinearAlgebra, LaTeXStrings, ProgressMeter, Statistics, Distributions
include("Lattice.jl")

printstyled("VEGETATION TURBIDITY LATTICE MODEL\n"; bold=true, underline=true, color=:light_blue)

# Define the lattice
n = 20
m = 20
N = n*m
L = Lattice(n, m)

# Initial condition of the lattice
x0 = 1.0::Float64.*ones(Float64, length(L.grid))
# Initial condition of the bifurcation parameter (degree of smooth muscle activation)
c = 3.0::Float64 # background turbidity
push!(x0, c)
# SDE's parameters
σ = 0.001::Float64 # noise level
push!(x0, σ)
ε = 0.01::Float64 # parameter's drift
push!(x0, ε)

# Model's parameters
R = 0.2::Float64 # dispersion rate
rv = 0.5::Float64 # maximum growth rate 
hv = 0.2::Float64 # half-saturation cover
r = rand(Uniform(0.6::Float64, 1.0::Float64), N) # half-saturation turbidity

# Sliding window's parameters
T = 10.00
δt = 1e-3
width = 20

# Definition of the SDE - deterministic part
function iip_det!(f, x, neigh, t)
        for j=1:N
                E = hv*x[N+1]/(hv + x[j])
                T = (r[j]^4 + E^4)/(r[j]^4) 
                f[j] = rv*x[j]*(1.0::Float64-x[j]*T) + R*(x[neigh[1,j]] + x[neigh[2,j]] + x[neigh[3,j]] + x[neigh[4,j]] - 4.0::Float64*x[j]) 
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
solution = solve(problem, EM(), dt=δt)
time = solution.t
states = solution.u

# Covariance analysis
EWS = zeros(Float64, (length(time)-width))
Mode = Array{Float64,2}(undef, N, (length(time)-width))
println("Performing the analysis of the covariance matrix modes")
@showprogress for t=1:(length(time)-width)
        # Get the states of the system across the specified time window 
        state = states[t:((t-1)+width)]
        # Assemble the snasphot matrix for the selected window
        ctr = 1
        X = Array{Float64, 2}(undef, width, N)
        for snapshot in state
                X[ctr, :] = snapshot[1:N]
                ctr +=1
        end
        # Compute the covariance matrix of the snapshots 
        Σ = cov(X, corrected=false)
        # Perform the eigendecomposition
        Λ = eigvals(Σ)
        Q = eigvecs(Σ)
        # Get leading eigenvalue and eigenvector as spatial EWS
        EWS[t] = (findmax(real(Λ)))[1]
        Mode[:,t] = real(Q[:,(findmax(real(Λ)))[2]])
end

# Postprocessing of the EWS using Python's timeseries classes 
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
        ts.detrend(mode='EMD', order=10)
        imfs = ts.imfs

        np.savetxt('./imfs1.csv', imfs[:,(imfs.shape)[1]-1], delimiter=',')
        np.savetxt('./imfs2.csv', imfs[:,(imfs.shape)[1]-2], delimiter=',')
        np.savetxt('./imfs3.csv', imfs[:,(imfs.shape)[1]-3], delimiter=',')
        np.savetxt('./imfs4.csv', imfs[:,(imfs.shape)[1]-4], delimiter=',')
        np.savetxt('./imfs5.csv', imfs[:,(imfs.shape)[1]-5], delimiter=',')
        np.savetxt('./imfs6.csv', imfs[:,(imfs.shape)[1]-6], delimiter=',')
        np.savetxt('./imfs7.csv', imfs[:,(imfs.shape)[1]-7], delimiter=',')

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


        plt.savefig("../results/lattice/vegetation_turbidity/fig2.9.png", transparent=True, dpi=600)
"""
# Analyse and detrend timeseries data
pyanalyse = py"analyse"
pyanalyse(EWS)

# Export images of the lattice time evolution
show(L, states, EWS, Mode)
