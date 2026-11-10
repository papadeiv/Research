"""
    Processing

Routines to postprocess and analyse the results of a simulation.
"""

# Estimate the return rate of a sample 
function estimate_return_rate(x)
  # Detrend the sample
  z = x .- mean(x)

  # Compute the variance and autocorrelation of the residuals
  γ = var(x)
  ρ = sum(z[1:(end-1)] .* z[2:end])/sum(abs2, z)

  # Compute the solution of the ordinary least square problem
  y = diff(x)./Δ
  θ = sum(z[1:(end-1)].*y)/sum(abs2, z[1:(end-1)])

  # Compute the estimators
  α_OUP = (σ^2/(2*γ))^2
  α_AC1 = ρ > 0 ? (log(ρ)/Δ)^2 : NaN
  α_OLS = θ^2

  return (OUP = α_OUP,
          AC1 = α_AC1,
          OLS = α_OLS)
end

# Compute the bias in the estimators above 
function compute_bias(μ)
  # Compute the autocorrelation at lag-1
  ρ = exp(λ(μ)*Δ)

  # Compute the bias in the leading eigenvalue
  b_OUP = -2*abs(λ(μ))*(1 + ρ + 2*ρ^2)/((Nt - 1)*(1 - ρ^2))
  b_AC1 = (-(5 + 2/ρ + 1/ρ^2)/(2*Δ*(Nt - 1)))
  b_OLS = Δ*λ(μ)^2/2 - (1 + 3*ρ)/(Δ*(Nt - 1))

  # Compute the uncertainty in the leading eigenvalue
  u = 2*abs(λ(μ))/(Δ*(Nt - 1))

  # Compute the bias in the return rate
  return (OUP = b_OUP^2 + 2*λ(μ)*b_OUP + u,
          AC1 = b_AC1^2 + 2*λ(μ)*b_AC1 + u,
          OLS = b_OLS^2 + 2*λ(μ)*b_OLS + u)
end

# Compute the mean, interquartile and extrema of the return rate estimators
function compute_summary_statistics(M)
  # Initialize empty data structures
  low_bound = fill(NaN, size(M, 1))
  first_quant = fill(NaN, size(M, 1))
  second_quant = fill(NaN, size(M, 1))
  third_quant = fill(NaN, size(M, 1))
  high_bound = fill(NaN, size(M, 1))

  for i in axes(M, 1)
    r = filter(isfinite, M[i,:])
    isempty(r) && continue
    low_bound[i] = quantile(r, 0.05)
    first_quant[i] = quantile(r, 0.25)
    second_quant[i] = mean(r)
    third_quant[i] = quantile(r, 0.75)
    high_bound[i] = quantile(r, 0.95)
  end

  return (low_bound = low_bound,
          first_quant = first_quant,
          second_quant = second_quant,
          third_quant = third_quant,
          high_bound = high_bound)
end
