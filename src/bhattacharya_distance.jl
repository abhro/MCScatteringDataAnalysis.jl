"""
    bcdistance(N₁, N₂)

Bhattacharya distance between two normal distributions
"""
function bcdistance(N₁, N₂)
    μ₁, σ₁ = params(N₁)
    μ₂, σ₂ = params(N₂)
    r_σ = (σ₁/σ₂)^2 # ratio of standard deviations
    # distance based on standard deviations, and means, respectively.
    distance_σ = log((r_σ + 1/r_σ + 2) / 4)
    distance_μ = (μ₁ - μ₂)^2/(σ₁^2 + σ₂^2)
    return 1/4 * (distance_σ + distance_μ)
end
