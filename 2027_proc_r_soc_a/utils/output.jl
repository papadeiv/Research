"""
Export data in csv format and saves figures.

Author: Davide Papapicco
Affil: U. of Auckland
Date: 02-12-2025
"""

function writeout(data, filename; path="../../res/data/")
        # Create the export directory if it does not exists
        fullpath = path * filename 
        mkpath(dirname(fullpath))

        # Write the input data to a csv file
        CSV.write(fullpath, Tables.table(data), delim=',', writeheader=false)
end

function savefig(path, fig)
        # Create the export directory if it doesn't exist
        fullpath = "../../res/fig/" * path 
        mkpath(dirname(fullpath))

        # Export the figure
        save(fullpath, fig)
end
