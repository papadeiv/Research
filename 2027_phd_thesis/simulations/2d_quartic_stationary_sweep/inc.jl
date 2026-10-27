using CairoMakie, Makie.Colors, LaTeXStrings
using ProgressMeter, Polynomials, Statistics
using DataFrames, CSV, MAT

# Avoid re-loading modules 
if !isdefined(Main, :DataInterface)
        include("../../src/DataInterface.jl")
        using .DataInterface
end

if !isdefined(Main, :SystemAnalysis)
        include("../../src/SystemAnalysis.jl")
        using .SystemAnalysis
end

if !isdefined(Main, :StatisticalMethods)
        include("../../src/StatisticalMethods.jl")
        using .StatisticalMethods
end

if !isdefined(Main, :PlottingTools)
        include("../../src/PlottingTools.jl")
        using .PlottingTools
end
