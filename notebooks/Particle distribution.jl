### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ a0342ed0-b9bc-11f0-aa1c-c52cb4db038a
using DrWatson

# ╔═╡ 0ca8a778-3e27-4350-9c6b-ec60c4bbdbb3
@quickactivate "MCScatteringDataAnalysis"

# ╔═╡ 39cf6f4b-e67c-472f-a00f-66fce4953228
using Unitful

# ╔═╡ f3a446b6-7fd4-4a0c-9551-c746f8a32c92
using Unitful: mp, c, yg

# ╔═╡ 9825cd92-8576-4554-90c0-97654c032526
using WGLMakie

# ╔═╡ 7fab5df1-c431-4aa1-83f2-55244a50175b
md"""
Particle distribution function from Warren et al. (2021)
```math
\frac{dN}{dp} = K p^{-σ} e^{-(p/p_\max)^χ}
```
"""

# ╔═╡ 784cdba7-4e74-4f69-9bc3-a33a8033c2f1
md"""
For brevity, define ``n(p) = dN/dp``.
"""

# ╔═╡ 5fe5684b-1a52-40ba-836c-a4cd35ea969f
md"""
A 'good' SI-ish unit for proton mass is yocto-grams
"""

# ╔═╡ 33dfa7f2-df32-4dc4-89ad-c45e7f20568a
u"mp" |> u"yg"

# ╔═╡ ecfe400a-1475-44ac-a1ba-9c0510fa3d94
n(p; σ, p_max, χ, p₀ = p_max, K = 1 / p₀) = K * (p / p₀)^-σ * exp(-(p / p_max)^χ)

# ╔═╡ c4a62cfc-9c4e-4e3c-b846-42f504f39526
ps = range(1.0e-4, 10^4, length = 300000) * u"mp*c" .|> u"yg*c"

# ╔═╡ 623454db-9f8e-48cd-b730-01086e1bf4d9
σ = 2.3

# ╔═╡ cd384e72-c781-46be-a7b1-b2b46c2ce2fc
p_max = 10^5 * u"mp*c" |> u"yg*c"

# ╔═╡ 30268d0c-ae65-4161-82e5-d8824dde1ee6
χ = 1.86;

# ╔═╡ 14c6c243-3a48-4071-bead-6fb8f6bf85df
ns = n.(ps; σ, p_max, χ)

# ╔═╡ 69b3e202-3d55-48a1-8181-c83e6b36c70e
let fig = Figure()
    ax = Axis(
        fig[1, 1],
        xlabel = "p",
        ylabel = "dN/dp",
        # dim2_conversion = Makie.UnitfulConversion(units_in_label=false)
    )

    lines!(ax, ps, ns)

    fig
end

# ╔═╡ 3f9ce62d-602d-4031-8bea-3803174febdd
let fig = Figure()
    ax = Axis(
        fig[1, 1],
        xlabel = "p",
        ylabel = "dN/dp",
        xscale = log10,
        yscale = log10,
        # dim1_conversion = Makie.UnitfulConversion(units_in_label=false),
        # dim2_conversion = Makie.UnitfulConversion(units_in_label=false),
    )

    lines!(ax, NoUnits.(ps ./ (mp * c)), NoUnits.(ns * (mp * c)))

    fig
end

# ╔═╡ b3a7a249-25fb-4f64-9d60-8aa7df20fe71
"""Maxwell-Juettner pdf, normalized against proton mass"""
f_MJ(p; θ) = p^2 * exp(-√(1 + (p / (mp * c))^2) / θ) / (mp * c)^3 |> unit(1 / p)

# ╔═╡ 21efacb8-0de4-4ff8-8e58-9705de029350
mp

# ╔═╡ 03b9fc61-8abc-4b37-8f5a-8e0b78cc743d
fs = f_MJ.(ps; θ = 30)

# ╔═╡ 6df48d91-72d5-4273-b3bb-ae66bd6f17b9
f_MB(v; T) = v^2 * exp(-v^2 / T)

# ╔═╡ 5cad4f59-5e92-496f-bb15-4366cc3c6d2b
lines(ustrip.(ps), f_MB.(ustrip.(ps); T = 90), axis = (; xscale = log10))

# ╔═╡ fd250922-abde-4c0a-8ef3-8f3c38a08072
let fig = Figure()
    ax = Axis(
        fig[1, 1],
        xlabel = "p",
        ylabel = "f(p)",
        xscale = log10,
        yscale = log10,
        xminorgridvisible = true,
        # dim1_conversion = Makie.UnitfulConversion(units_in_label=false),
        # dim2_conversion = Makie.UnitfulConversion(units_in_label=false),
    )

    lines!(ax, ustrip.(ps), ustrip.(fs))

    fig
end

# ╔═╡ Cell order:
# ╠═a0342ed0-b9bc-11f0-aa1c-c52cb4db038a
# ╠═0ca8a778-3e27-4350-9c6b-ec60c4bbdbb3
# ╠═39cf6f4b-e67c-472f-a00f-66fce4953228
# ╠═f3a446b6-7fd4-4a0c-9551-c746f8a32c92
# ╟─7fab5df1-c431-4aa1-83f2-55244a50175b
# ╟─784cdba7-4e74-4f69-9bc3-a33a8033c2f1
# ╟─5fe5684b-1a52-40ba-836c-a4cd35ea969f
# ╠═33dfa7f2-df32-4dc4-89ad-c45e7f20568a
# ╠═9825cd92-8576-4554-90c0-97654c032526
# ╠═ecfe400a-1475-44ac-a1ba-9c0510fa3d94
# ╠═c4a62cfc-9c4e-4e3c-b846-42f504f39526
# ╠═14c6c243-3a48-4071-bead-6fb8f6bf85df
# ╠═623454db-9f8e-48cd-b730-01086e1bf4d9
# ╠═cd384e72-c781-46be-a7b1-b2b46c2ce2fc
# ╠═30268d0c-ae65-4161-82e5-d8824dde1ee6
# ╠═69b3e202-3d55-48a1-8181-c83e6b36c70e
# ╠═3f9ce62d-602d-4031-8bea-3803174febdd
# ╠═b3a7a249-25fb-4f64-9d60-8aa7df20fe71
# ╠═21efacb8-0de4-4ff8-8e58-9705de029350
# ╠═03b9fc61-8abc-4b37-8f5a-8e0b78cc743d
# ╠═6df48d91-72d5-4273-b3bb-ae66bd6f17b9
# ╠═5cad4f59-5e92-496f-bb15-4366cc3c6d2b
# ╠═fd250922-abde-4c0a-8ef3-8f3c38a08072
