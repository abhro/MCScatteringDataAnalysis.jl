"""
    bcdistance(N₁, N₂)

Bhattacharya distance between two normal distributions
"""
function bcdistance(N₁, N₂)
    μ₁, σ₁ = params(N₁)
    μ₂, σ₂ = params(N₂)
    r_σ = (σ₁ / σ₂)^2 # ratio of standard deviations
    # distance based on standard deviations, and means, respectively.
    distance_σ = log((r_σ + 1/r_σ + 2) / 4)
    distance_μ = (μ₁ - μ₂)^2 / (σ₁^2 + σ₂^2)
    return (distance_σ + distance_μ) / 4
end

"""
    bcdistances(NV₁, NV₂)

Given two vectors of distributions, return the Bhattacharya distance between
each pair of distributions in the two vectors. (May need to rewrite this
sentence for clarity.)
If either of the two distributions is missing, then that distance will instead
contain missing.
"""
function bcdistances(NV₁, NV₂)
    n = length(NV₁)
    if n != length(NV₂)
        throw(DimensionMismatch("The two vectors of distributions must have the same length"))
    end
    distances = zeros(Union{Missing, Float64}, n)
    for (idx, (distrib1, distrib2)) in enumerate(zip(NV₁, NV₂))
        if ismissing(distrib1) || ismissing(distrib2)
            distances[idx] = missing
        else
            distances[idx] = bcdistance(distrib1, distrib2)
        end
    end
    return distances
end
