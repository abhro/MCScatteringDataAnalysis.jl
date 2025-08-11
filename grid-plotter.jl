### A Pluto.jl notebook ###
# v0.20.4

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
## Import packages and set up base variables
"""

# ╔═╡ cd50f95b-8e0b-44a1-982b-52f08192d7b4
TableOfContents(depth = 6)

# ╔═╡ 95f0d053-cf1f-41e6-8a60-437ad4b547f1


# ╔═╡ 41051c07-b2bc-4890-ad4e-9261f2f6c7ed
md"""
Increase cell width
"""

# ╔═╡ 846ce513-17b4-49aa-8823-15e3f90c7dcb
html"""<style>
main {
    max-width: 70%;
    padding-left: max(360px, 10%);
    padding-right: 0%;
}
</style>"""

# ╔═╡ 9bd001d5-d346-4f92-890c-358fce5e64ea
md"""
# Read data files
"""

# ╔═╡ ee889b66-662c-4e86-a42e-50e02c99323b
md"""
## Create data file filter
"""

# ╔═╡ 6eefa524-9d4a-418e-b0fd-9a770986723c
grid_df = load_object("grid.jld2")

# ╔═╡ 6e28bd19-e152-420b-905b-58ffcfc84f8f
describe(grid_df)

# ╔═╡ 3213f4b0-a5f7-4fa3-a53b-c45a620953e7


# ╔═╡ 7205a1dd-3e4b-40d4-9584-e10607a6139f
md"""
# Constants and functions
"""

# ╔═╡ dd93decc-706c-4c96-8e60-2ac734916e81


# ╔═╡ Cell order:
# ╠═8cecd3b0-d2b1-11ef-35c7-1fb9d13a3262
# ╟─8cecd3b0-d2b1-11ef-2fa6-932b3621a372
# ╟─50eff28b-b1ba-4f29-94d3-0d8862952cd1
# ╠═a0c6b97e-f294-4c6e-ac80-7bfc902d58e7
# ╠═87c538d7-28ab-40ef-84d7-1b34feceec44
# ╠═e4034545-f28d-442e-9e45-9d57d6a07ef3
# ╠═cd50f95b-8e0b-44a1-982b-52f08192d7b4
# ╠═95f0d053-cf1f-41e6-8a60-437ad4b547f1
# ╟─41051c07-b2bc-4890-ad4e-9261f2f6c7ed
# ╟─846ce513-17b4-49aa-8823-15e3f90c7dcb
# ╟─9bd001d5-d346-4f92-890c-358fce5e64ea
# ╟─ee889b66-662c-4e86-a42e-50e02c99323b
# ╠═6eefa524-9d4a-418e-b0fd-9a770986723c
# ╠═6e28bd19-e152-420b-905b-58ffcfc84f8f
# ╠═3213f4b0-a5f7-4fa3-a53b-c45a620953e7
# ╟─7205a1dd-3e4b-40d4-9584-e10607a6139f
# ╠═dd93decc-706c-4c96-8e60-2ac734916e81
