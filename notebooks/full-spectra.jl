### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ c499612e-379f-11f1-b64f-393d1cf53318
import Pkg; Pkg.activate(Base.current_project())

# ╔═╡ 43e50caf-9893-4e0c-88d8-09abf61e2d64
using CairoMakie; set_theme!(theme_latexfonts())

# ╔═╡ c499612e-379f-11f1-ab84-331adab584d1
χ = 1.86;

# ╔═╡ 9a9e7992-283e-4f89-ab4d-539344068c9a
p_c = 1;    # where we switch from thermal to power law

# ╔═╡ e7042973-785c-4649-8368-0569778d1fe0
p_max = 1800; # exponential decay scale

# ╔═╡ 790b8d61-dd03-4f69-b7ef-01fa1453cfa3
μ = 50;     # relativistic coldness, mc²/kT

# ╔═╡ 0ea7b513-4e73-4bed-8955-57fe37c65009
σ = 4.0;    # power-law spectral index

# ╔═╡ 224d49a5-d41a-4c3b-8d6f-df12f61227ab
function distribution(p; p_c, p_max, χ, μ, σ)
    C_TH = p_c^(-(σ + 2)) * exp(-(p_c / p_max)^χ + μ * sqrt(1 + p_c^2))
    if p ≤ p_c
        nₚ = C_TH * p^2 * exp(-μ * sqrt(1 + p^2))
    else
        nₚ = p^-σ * exp(-(p / p_max)^χ)
    end
    return nₚ
end

# ╔═╡ 6acfb0c8-b7ee-41d3-8a68-8bf42e5b916a
function log_distribution(p; p_c, p_max, χ, μ, σ)
    G_TH = -(σ + 2) * log10(p_c) - (p_c / p_max)^χ + μ * sqrt(1 + p_c^2)
    logp = log10(p)
    if p ≤ p_c
        lognₚ = G_TH + 2 * logp + -μ * sqrt(1 + p^2)
    else
        lognₚ = -σ * p + -(p / p_max)^χ
    end
    return lognₚ
end

# ╔═╡ 285bd427-4056-4f63-a166-d4435eca2455
p = logrange(0.01, 9000, length = 1000);

# ╔═╡ bfaeb834-9f55-4006-b84e-e2a78bd306da
nₚ = distribution.(p; p_c, p_max, χ, μ, σ) * 1.0e9;

# ╔═╡ 91c02096-25c4-401c-b305-e0c0cd6be583
let fig = Figure()
    ax = Axis(
        fig[1, 1];
        xscale = log10, yscale = log10,
        xlabel = L"p", ylabel = L"n_p",
        xminorticksvisible = true,
        yminorticksvisible = true,
        xminorgridvisible = true,
        yminorgridvisible = true,
    )
    lines!(ax, p, nₚ)
    fig
end

# ╔═╡ Cell order:
# ╠═c499612e-379f-11f1-b64f-393d1cf53318
# ╠═c499612e-379f-11f1-ab84-331adab584d1
# ╠═43e50caf-9893-4e0c-88d8-09abf61e2d64
# ╠═9a9e7992-283e-4f89-ab4d-539344068c9a
# ╠═e7042973-785c-4649-8368-0569778d1fe0
# ╠═790b8d61-dd03-4f69-b7ef-01fa1453cfa3
# ╠═0ea7b513-4e73-4bed-8955-57fe37c65009
# ╠═224d49a5-d41a-4c3b-8d6f-df12f61227ab
# ╠═6acfb0c8-b7ee-41d3-8a68-8bf42e5b916a
# ╠═285bd427-4056-4f63-a166-d4435eca2455
# ╠═bfaeb834-9f55-4006-b84e-e2a78bd306da
# ╠═91c02096-25c4-401c-b305-e0c0cd6be583
