const log2π = log(2π)

function aic{F<:Real,R<:Real}(y::Vector{F}, y_hat::Vector{R}, p::Int)
  const n = length(y)
  n == length(y_hat) || throw(DimensionMismatch("y and y_hat have to be of the same size"))
  const r = y_hat - y
  n * (log2π + 1 + log(dot(r,r)/n)) + 2(p+1)
end

