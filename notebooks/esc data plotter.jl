### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ cd76c040-d2b5-11ef-29c3-8d8f9cb45061
using DrWatson

# ╔═╡ e5e0e4e2-2df1-4536-9cc5-bdcec6fc13de
@quickactivate "MCScatteringDataAnalysis"

# ╔═╡ b8029e47-60db-4574-af33-51d74b25bb60
using CSV, DataFrames

# ╔═╡ 0c4a7349-4ecb-454c-9c8d-4dcf2b66ba94
using PlutoUI

# ╔═╡ 590b9a88-33ff-4400-82ba-40be75d8da23
using CairoMakie

# ╔═╡ cd76c040-d2b5-11ef-3898-079cb259eed1
md"""
# Preamble
"""

# ╔═╡ 9e48bcdc-b9d3-4a71-81bf-d45a3a1f4037
md"""
## Import packages
"""

# ╔═╡ eb9c48c8-42a5-4416-a1c0-ff43c7e91aae
md"""
## Configure notebook appearance
"""

# ╔═╡ e950378a-3c65-4148-ad9d-7f1f82afaa01
TableOfContents(depth = 6)

# ╔═╡ 194f6821-15d4-4398-8daa-81164ec55578
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

# ╔═╡ ddefdc26-3cdc-480d-a5d5-32dc7c39d4a1
md"""
# Read data files
"""

# ╔═╡ 8f8beeee-3b47-4463-a233-1024478a8cf8
escdf = CSV.read(datadir("Lorentz-5-processed", "dNdp-esc.csv.gz"), DataFrame);

# ╔═╡ 61632ddf-55a6-4fd5-84c0-8344b007c70d
escdf

# ╔═╡ f802e55a-d945-4c88-a285-8d3b906211f2
md"""
# Plot data
"""

# ╔═╡ da1ab6df-120a-4afd-a532-70faaa10dcc4
let f = Figure()
    ax1 = Axis(f[1,1], title = "dN/dp, escaping, upstream")
    ax2 = Axis(f[2,1], title = "dN/dp, escaping, downstream", xlabel = "log(p)")

    for df in groupby(escdf, [:initial_seed, :itr])
        scatter!(ax1, df.log_p_cgs, df.log_dNdp_esc_UpS_IF, label = "Upstream")
        scatter!(ax2, df.log_p_cgs, df.log_dNdp_esc_DwS_IF, label = "Downstream")
    end

    axislegend(ax1)
    axislegend(ax2)

    f
end

# ╔═╡ Cell order:
# ╟─cd76c040-d2b5-11ef-3898-079cb259eed1
# ╟─9e48bcdc-b9d3-4a71-81bf-d45a3a1f4037
# ╠═cd76c040-d2b5-11ef-29c3-8d8f9cb45061
# ╠═e5e0e4e2-2df1-4536-9cc5-bdcec6fc13de
# ╠═b8029e47-60db-4574-af33-51d74b25bb60
# ╠═0c4a7349-4ecb-454c-9c8d-4dcf2b66ba94
# ╠═590b9a88-33ff-4400-82ba-40be75d8da23
# ╟─eb9c48c8-42a5-4416-a1c0-ff43c7e91aae
# ╠═e950378a-3c65-4148-ad9d-7f1f82afaa01
# ╟─194f6821-15d4-4398-8daa-81164ec55578
# ╟─ddefdc26-3cdc-480d-a5d5-32dc7c39d4a1
# ╠═8f8beeee-3b47-4463-a233-1024478a8cf8
# ╠═61632ddf-55a6-4fd5-84c0-8344b007c70d
# ╟─f802e55a-d945-4c88-a285-8d3b906211f2
# ╠═da1ab6df-120a-4afd-a532-70faaa10dcc4
