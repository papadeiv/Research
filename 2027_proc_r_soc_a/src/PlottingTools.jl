module PlottingTools

# Import packages
using CairoMakie, Makie.Colors, LaTeXStrings
using DocStringExtensions

# Increase the definition of the figures
CairoMakie.activate!(; px_per_unit = 2)

# Export namespaces
export savefig 

end # module
