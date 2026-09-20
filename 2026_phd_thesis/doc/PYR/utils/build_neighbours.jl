function build_neighbours(I::Tuple{Int,Int})
    neigh = Array{Int,2}(undef,4,prod(I))
    rows = I[1]
    cols = I[2]
    ctr = 0
    for ix in 1:cols
        for iy in 1:rows
            ctr += 1
            N = mod1(iy - 1, rows), ix 
            E = iy, mod1(ix + 1, cols)
            S = mod1(iy + 1, rows), ix
            W = iy, mod1(ix - 1, cols)
            neigh[1,ctr] = N[1] + I[1]*(N[2]-1)
            neigh[2,ctr] = W[1] + I[1]*(W[2]-1)
            neigh[3,ctr] = S[1] + I[1]*(S[2]-1)
            neigh[4,ctr] = E[1] + I[1]*(E[2]-1)
        end
    end
    return neigh
end
