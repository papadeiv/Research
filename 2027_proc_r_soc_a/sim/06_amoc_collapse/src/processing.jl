"""
    Processing

Routines to estimate the return rate and extrapolate the tipping time.
"""

# Linear detrending of a non-stationary timeseries
function detrend(t, x)
  B = hcat(ones(length(t)), t .- mean(t))
  return x .- B*(B\x)
end

# Estimate the return rate
function estimate_return_rate(t, x)
  # Locate the first and last observation in the admissible time horizon
  T1_idx = searchsortedfirst(t, T1)
  T2_idx = searchsortedlast(t, T2)

  # Initialize empty arrays
  timesteps = Float64[]
  rate = Float64[]
  uncertainty = Float64[]

  # Loop over the windows
  for idx in T1_idx:Nw:(T2_idx - Nw + 1)
    # Obtain windowed residuals 
    window = idx:(idx + Nw - 1)
    z = detrend(t[window], x[window])

    # Estimate the leading eigenvalue and its uncertainty 
    y = diff(z)./Δ
    λ = sum(z[1:(end-1)].*y)/sum(abs2, z[1:(end-1)])
    u = 2*abs(λ)/(Δ*(Nw - 1))

    # Compute the return rate and its uncertainty 
    α = λ^2 - u
    u = 4*(λ^2)*u

    # Store the estimate at the midpoint of the window
    push!(timesteps, t[idx + Nw ÷ 2])
    push!(rate, α)
    push!(uncertainty, u)
  end

  return (t = timesteps, α = rate, u = uncertainty)
end

# Extrapolate the breaking time
function extrapolate_breaking_time(estimates)
  # Filter NaN and Inf 
  valid = isfinite.(estimates.α) .& isfinite.(estimates.u) .& (estimates.u .> 0)
  count(valid) < 2 && return (τ = NaN, q = NaN, m = NaN)
  t = estimates.t[valid]
  α = estimates.α[valid]
  u = estimates.u[valid]

  # Assemble the model and the weight matrix
  B = hcat(ones(length(t)), t)
  W = Diagonal(1.0./u)

  # Solve the weighted normal equation
  q, m = (B'*W*B) \ (B'*W*α)

  return (τ = -q/m, q = q, m = m)
end
