### A Pluto.jl notebook ###
# v0.20.15

using Markdown
using InteractiveUtils

# ╔═╡ 8cecd3b0-d2b1-11ef-35c7-1fb9d13a3262
import Pkg; Pkg.activate(Base.current_project())

# ╔═╡ a0c6b97e-f294-4c6e-ac80-7bfc902d58e7
using JLD2, DataFrames

# ╔═╡ 87c538d7-28ab-40ef-84d7-1b34feceec44
using PlutoUI

# ╔═╡ e4034545-f28d-442e-9e45-9d57d6a07ef3
using CairoMakie

# ╔═╡ 8cecd3b0-d2b1-11ef-2fa6-932b3621a372
md"""
# Preamble
"""

# ╔═╡ 50eff28b-b1ba-4f29-94d3-0d8862952cd1
md"""
## Import packages
"""

# ╔═╡ f6b1dddc-6cfc-4cac-b7ac-4c476205793f
md"""
## Configure notebook appearance
"""

# ╔═╡ cd50f95b-8e0b-44a1-982b-52f08192d7b4
TableOfContents(depth = 6)

# ╔═╡ 846ce513-17b4-49aa-8823-15e3f90c7dcb
# Increase cell width
html"""
<style>
main {
    max-width: 83%;
    padding-left: max(300px, 5%);
    padding-right: 0%;
}
</style>
"""

# ╔═╡ 9bd001d5-d346-4f92-890c-358fce5e64ea
md"""
# Read data files
"""

# ╔═╡ 5deeedc5-7c11-495f-a3db-27f2e6f5dfd9
const datadir = "G:/My Drive/MC Scattering/Processed-data/";

# ╔═╡ 6eefa524-9d4a-418e-b0fd-9a770986723c
grid_df = load_object(joinpath(datadir, "grid.jld2"))

# ╔═╡ 6e28bd19-e152-420b-905b-58ffcfc84f8f
describe(grid_df)

# ╔═╡ Cell order:
# ╟─8cecd3b0-d2b1-11ef-2fa6-932b3621a372
# ╟─50eff28b-b1ba-4f29-94d3-0d8862952cd1
# ╠═8cecd3b0-d2b1-11ef-35c7-1fb9d13a3262
# ╠═a0c6b97e-f294-4c6e-ac80-7bfc902d58e7
# ╠═87c538d7-28ab-40ef-84d7-1b34feceec44
# ╠═e4034545-f28d-442e-9e45-9d57d6a07ef3
# ╟─f6b1dddc-6cfc-4cac-b7ac-4c476205793f
# ╠═cd50f95b-8e0b-44a1-982b-52f08192d7b4
# ╟─846ce513-17b4-49aa-8823-15e3f90c7dcb
# ╟─9bd001d5-d346-4f92-890c-358fce5e64ea
# ╠═5deeedc5-7c11-495f-a3db-27f2e6f5dfd9
# ╠═6eefa524-9d4a-418e-b0fd-9a770986723c
# ╠═6e28bd19-e152-420b-905b-58ffcfc84f8f
