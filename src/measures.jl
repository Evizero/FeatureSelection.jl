const log2π = log(2π)

"""
The Akaike information criterion
"""
function aic{F<:Real,R<:Real}(y::Vector{F}, ŷ::Vector{R}, p::Int, k::FloatingPoint=2.)
  const n = length(y)
  n == length(ŷ) || throw(DimensionMismatch("y and y_hat have to be of the same size"))
  const r = ŷ - y
  n * (log2π + 1 + log(dot(r,r)/n)) + k*(p+1)
end

"""
The Bayesian information criterion, also known as Schwarz's Bayesian criterion.
"""
function bic{F<:Real,R<:Real}(y::Vector{F}, ŷ::Vector{R}, p::Int)
  const n = length(y)
  n == length(ŷ) || throw(DimensionMismatch("y and y_hat have to be of the same size"))
  const k = log(n)
  const r = ŷ - y
  n * (log2π + 1 + log(dot(r,r)) - k) + k*(p+1)
end

