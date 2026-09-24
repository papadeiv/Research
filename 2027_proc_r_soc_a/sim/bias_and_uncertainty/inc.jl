using CairoMakie, Makie.Colors, LaTeXStrings
using LinearAlgebra, Statistics  
using ProgressMeter

# Avoid re-loading DataInterface 
if !isdefined(Main, :DataInterface)
        include("../../src/DataInterface.jl")
        using .DataInterface
end

# Avoid re-loading SystemAnalysis 
if !isdefined(Main, :SystemAnalysis)
        include("../../src/SystemAnalysis.jl")
        using .SystemAnalysis
end

# Avoid re-loading StatisticalMethods
if !isdefined(Main, :StatisticalMethods)
        include("../../src/StatisticalMethods.jl")
        using .StatisticalMethods
end

# Avoid re-loading PlottingTools
if !isdefined(Main, :PlottingTools)
        include("../../src/PlottingTools.jl")
        using .PlottingTools
end
