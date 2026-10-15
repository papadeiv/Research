"""
    Processing

Routines to postprocess and analyse the results of a simulation.
"""

# Detrend a window by removing its linear trend
function detrend(t, x)
  B = hcat(ones(length(t)), t .- mean(t))
  return x .- B*(B\x)
end

# Estimate the leading eigenvalue using (3.7), (3.11) and (3.16)
function estimate_eigenvalue(z)
  # Compute the variance and autocorrelation of the residuals
  γ = var(z)
  ρ = sum(z[1:(end-1)] .* z[2:end])/sum(abs2, z)

  # Compute the solution of the ordinary least square problem
  y = diff(z)./Δ
  θ = sum(z[1:(end-1)].*y)/sum(abs2, z[1:(end-1)])

  return (OUP = -σ^2/(2*γ),
          AC1 = ρ > 0 ? log(ρ)/Δ : NaN,
          OLS = θ)
end

# Correct the return rate for the uncertainty in the eigenvalue using (4.13) and (4.14)
function estimate_return_rate(λ̂)
  # Compute the uncertainty in the leading eigenvalue
  uλ = 2*abs(λ̂)/(Δ*(Nw - 1))

  # Compute the corrected return rate and its uncertainty
  return (α = λ̂^2 - uλ,
          u = 4*λ̂^2*uλ)
end

# Extrapolate the breaking time by weighted least squares
function estimate_breaking_time(t, α, u)
  # Discard the windows without an estimate
  k = findall(i -> isfinite(α[i]) && isfinite(u[i]) && u[i] > 0, eachindex(α))
  length(k) < 2 && return (τ = NaN, q = NaN, m = NaN)

  # Assemble the model and the weight matrix
  B = hcat(ones(length(k)), t[k])
  W = Diagonal(1.0./u[k])

  # Solve the weighted normal equations
  q, m = (B'*W*B) \ (B'*W*α[k])

  return (τ = -q/m, q = q, m = m)
end

# Compute the ensemble median of the estimates in every window
function compute_median(M)
  rows = [filter(isfinite, M[i,:]) for i in axes(M, 1)]
  return [isempty(r) ? NaN : median(r) for r in rows]
end

# Summarise a sample by its thinned path and the return rates estimated in every window
function summarise_sample(path)
  # Initialize empty data structures
  estimates = (OUP = fill(NaN, length(tw)),
               AC1 = fill(NaN, length(tw)),
               OLS = fill(NaN, length(tw)))
  uncertainties = (OUP = fill(NaN, length(tw)),
                   AC1 = fill(NaN, length(tw)),
                   OLS = fill(NaN, length(tw)))

  # Loop over the rolling windows
  for (i, s) in enumerate(starts)
    # Skip the windows beyond the escape of the sample
    s + Nw - 1 <= length(path.u) || continue

    # Detrend the sample within the window
    w = s:(s + Nw - 1)
    z = detrend(path.t[w], path.u[w])

    # Estimate the leading eigenvalue and the return rate of the window
    λ̂ = estimate_eigenvalue(z)
    OUP = estimate_return_rate(λ̂.OUP)
    AC1 = estimate_return_rate(λ̂.AC1)
    OLS = estimate_return_rate(λ̂.OLS)
    estimates.OUP[i], uncertainties.OUP[i] = OUP.α, OUP.u
    estimates.AC1[i], uncertainties.AC1[i] = AC1.α, AC1.u
    estimates.OLS[i], uncertainties.OLS[i] = OLS.α, OLS.u
  end

  return (t = thin(path.t, Nt), x = thin(path.u, Nt),
          estimates = estimates, uncertainties = uncertainties)
end

# Thin the ensemble down to the solutions that are drawn
function collect_solutions(Ω)
  return [(t = sample.t, x = sample.x) for sample in Ω[1:max(1, Ne ÷ Np):end]]
end

# Print the median and interquartile range of the breaking times
function print_breaking_times(breaking_times)
  @printf("τ = %.1f, τpb = %.1f, %d windows of %d samples\n", τ, τpb, length(tw), Nw)
  for key in (:OUP, :AC1, :OLS)
    τ̂ = filter(isfinite, getproperty(breaking_times, key))
    @printf("  %s: median %.1f (%+.2f%%), IQR [%.1f, %.1f]\n", key, median(τ̂),
            100*(median(τ̂)/τ - 1), quantile(τ̂, 0.25), quantile(τ̂, 0.75))
  end
end

# Compute the plotting range of the breaking times between the 1st and 99th percentiles, widened to include τ
function compute_range(τ̂)
  low, high = quantile(τ̂, 0.01), quantile(τ̂, 0.99)
  pad = 0.05*(high - low)

  return (min(low - pad, τ - pad), max(high + pad, τ + pad))
end

# Thin a vector down to at most n points
thin(v, n) = v[round.(Int, range(1, length(v), length=min(n, length(v))))]
