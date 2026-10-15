"""
    Processing

Routines to postprocess and analyse the results of a simulation.
"""

# Detrend a window by removing its linear trend
function detrend(t, x)
  B = hcat(ones(length(t)), t .- mean(t))
  return x .- B*(B\x)
end

# Estimate the return rate using (3.16), corrected for the uncertainty using (4.13) and (4.14)
function estimate_return_rate(z)
  # Compute the solution of the ordinary least square problem
  y = diff(z)./Δ
  λ̂ = sum(z[1:(end-1)].*y)/sum(abs2, z[1:(end-1)])

  # Compute the uncertainty in the leading eigenvalue
  uλ = 2*abs(λ̂)/(Δ*(Nw - 1))

  return (α = λ̂^2 - uλ, u = 4*λ̂^2*uλ)
end

# Summarise a sample by its thinned path and the return rates estimated in every window
function summarise_sample(path)
  α, u = fill(NaN, length(tw)), fill(NaN, length(tw))

  # Loop over the rolling windows that the sample reaches
  for (i, s) in enumerate(starts)
    s + Nw - 1 <= length(path.u) || continue
    w = s:(s + Nw - 1)
    α[i], u[i] = estimate_return_rate(detrend(path.t[w], path.u[w]))
  end

  # Record the tipping time of the sample, if it tipped
  tip = path.t[end] < T - Δ ? path.t[end] : NaN

  return (t = thin(path.t, Nt), x = thin(path.u, Nt), α = α, u = u, tip = tip)
end

# Extrapolate the breaking time by weighted least squares
function estimate_breaking_time(t, α, u)
  # Discard the windows without an estimate
  k = findall(i -> isfinite(α[i]) && u[i] > 0, eachindex(α))
  length(k) < 2 && return (τ = NaN, q = NaN, m = NaN)

  # Solve the weighted normal equations
  B = hcat(ones(length(k)), t[k])
  W = Diagonal(1.0./u[k])
  q, m = (B'*W*B) \ (B'*W*α[k])

  return (τ = -q/m, q = q, m = m)
end

# Compute the median of the estimates in every window
compute_median(M) = [median(filter(isfinite, M[i,:])) for i in axes(M, 1)]

# Split a sample into the parts lying outside and inside the admissible horizon
function split_admissible(t)
  return (before = t .<= tw[1],
          admissible = tw[1] .<= t .<= tw[end],
          after = t .>= tw[end])
end

# Compute the tail integral of the fourth power of Ai by the trapezoidal rule
function airy_tail(w; wmax=10.0, n=3000)
  u = range(w, wmax, length=n)
  y = airyai.(u).^4

  return (sum(y) - (y[1] + y[end])/2)*step(u)
end

# Compute the small noise law of the tipping time, supported between the first two zeros of Bi
function p_small_noise(t)
  # Rescale the time and reject the points outside the support
  w = ε^(1/3)*(τ - t)
  -3.2711 < w < -1.1737 || return 0.0

  # Evaluate the Airy functions and their tail integral
  A, B, g = airyai(w), airybi(w), airy_tail(w)

  return ε^(5/6)/(sqrt(4π^5*D)*B^2*sqrt(g))*exp(-ε*A^2/(4π^2*D*B^2*g))
end

# Compute the number of histogram bins by Scott's rule
scott(x) = max(1, ceil(Int, (maximum(x) - minimum(x))/(3.49*std(x)*length(x)^(-1/3))))

# Thin a vector down to at most n points
thin(v, n) = v[round.(Int, range(1, length(v), length=min(n, length(v))))]
