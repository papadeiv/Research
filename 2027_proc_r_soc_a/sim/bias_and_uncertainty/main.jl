"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/figs.jl")

# Main algorithm 
function main()
        # Solve the ensemble problems and export the results
        tipped = generate_samples()

        # Check whether any sample path has tipped
        if tipped
                @goto stop
        else
                # Perform the statistical analysis for the current simulation setup
                analysis()
        end

        # Skip to the end
        @goto skip

        # Stop the analysis 
        @label stop
        printstyled("A sample path has tipped\n"; bold=true, underline=true, color=:light_blue)

        # End the simulation
        @label skip
end

# Execute the main 
main()
