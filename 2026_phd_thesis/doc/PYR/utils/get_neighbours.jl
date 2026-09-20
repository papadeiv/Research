function get_neighbours(Lattice::Lattice, node_index::Int64)
        neighbours = Lattice.connectivity[:,node_index]
        return neighbours
end
