"""
Read and format data from txt, csv and mat files.

Author: Davide Papapicco
Affil: U. of Auckland
Date: 02-12-2025
"""

function readin(filename; path="../../res/data/")
        # Give full path to the reader
        fullpath = path * filename 

        # Read the data in dataframe form
        df = DataFrame(CSV.File(fullpath; delim=',', header=false))

        # Define the time of datatype to be returned based on the number of columns in the dataframe
        if ncol(df) == 1 
                # Convert the dataframe in vector form
                global data = Vector{Float64}(undef, nrow(df))
                data = df[:,1]
        else
                # Convert the dataframe in matrix form
                data = Matrix{Float64}(undef, nrow(df), ncol(df))
                for n in 1:ncol(df)
                        data[:,n] = df[:,n]
                end
        end

        # Export the array
        return data 
end
