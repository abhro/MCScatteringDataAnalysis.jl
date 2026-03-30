### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 67d360be-2973-11f1-85b4-13446ccb6184
# ╠═╡ skip_as_script = true
#=╠═╡
import Pkg; Pkg.activate(Base.current_project())
  ╠═╡ =#

# ╔═╡ 67d360be-2973-11f1-81b0-e9f876bff0c9
using CairoMakie

# ╔═╡ e2eddf74-c76a-43e0-8c44-d1b06894c3cc
using Distributions

# ╔═╡ 2d19b46a-7de8-4156-b47b-1292642b5a54
using Random

# ╔═╡ 7d281913-c40d-4452-941f-7f5a28d62328
rng = Xoshiro(5)

# ╔═╡ df4b7c63-cb3e-4432-a4db-78f14f093a06
# compression ratio
r = 4;

# ╔═╡ 6030ac76-e906-46fa-a68f-20a41c7ff130
# number of upstream particles
n_upstream = 50;

# ╔═╡ 67d360be-2973-11f1-9cd7-4d87cb3d3c08
# number of downstream_particles
n_downstream = n_upstream * r

# ╔═╡ 8cf4fc55-daac-4915-ab30-6150fd948efd
y_dist = Uniform(-2, 2)

# ╔═╡ f2ac1dee-3b4a-4f75-b39e-5e29c457a9b2
x_upstream_dist = Uniform(-5, -0.1)

# ╔═╡ 808a2a5f-a6e3-4660-9118-c5a430d20c34
x_downstream_dist = Uniform(0.1, 5)

# ╔═╡ ce51675e-210d-4626-99fd-febf01108d67
begin
    # x and y coordinates for upstream particles
    x_upstream = rand(rng, x_upstream_dist, n_upstream)
    y_upstream = rand(rng, y_dist, n_upstream)

    # x and y coordinates for downstream particles
    x_downstream = rand(rng, x_downstream_dist, n_downstream)
    y_downstream = rand(rng, y_dist, n_downstream)
end;

# ╔═╡ 15823133-b986-4d03-be43-1d01a4aa74da
begin
    points_downstream = Point2.(x_downstream, y_downstream)
    points_upstream = Point2.(x_upstream, y_upstream)
end;

# ╔═╡ 041d4f92-c2e5-435a-a07f-31339554583e
begin
    # downstream particles that get a velocity vector
    v_downstream = sample(rng, points_downstream, 40)
    # upstream particles that get a velocity vector
    v_upstream = sample(rng, points_upstream, 30)
end;

# ╔═╡ dd966d52-0d58-4ce2-a1d2-a3d2d0e45662
# size of arrow for upstream velocity
u_upstream = 0.45

# ╔═╡ 1245e54f-c85c-445d-956c-117a01115126
# size of arrow for downstream velocity
u_downstream = u_upstream / √√r

# ╔═╡ 928a7e9c-9010-4eb4-87c9-10213f7a3164
markersize = 5;

# ╔═╡ 218cf991-e152-4b45-ae84-ee2d5d686e28
shaftwidth = 1;

# ╔═╡ 1b47fa68-d9d1-486d-88dd-2c094248b974
set_theme!(theme_latexfonts())

# ╔═╡ c6a99716-b7f9-4586-b0e8-257ded2564de
upstream_color = Makie.wong_colors()[1]

# ╔═╡ 36f34b4f-a3f3-4637-9db6-e6d435e54d18
downstream_color = Makie.wong_colors()[3]

# ╔═╡ bbc09185-bb2b-4bf0-84fc-3d644851c3f3
schematic_figure = let
    fig = Figure()
    ax = Axis(
        fig[1,1];
        xgridvisible = false, ygridvisible = false,
        xticks = -5:5,
        xticksvisible = true, xminorticksvisible = true,
        yticksvisible = false,
        # xticklabelsvisible = false,
        yticklabelsvisible = false,
        xlabel = L"x",
        # title = "Schematic of plasma flow in shock frame",
    )
    hidespines!(ax, :l, :t, :r)

    text!(ax, -4.6, 1.3, text = L"\textbf{upstream}, \; u_\text{u,sf}, \; n_\text{u,sf}",  color = upstream_color)
    text!(ax, 1, 2, text = L"\textbf{downstream}, \; u_\text{d,sf}, \; n_\text{d,sf}", color = downstream_color)

    arrows2d!(
        ax, v_upstream, fill(Point2(u_upstream, 0), length(v_upstream));
        color = upstream_color,
        # label = L"u_\text{u,sf}",
        shaftwidth,
        tiplength = 4, tipwidth = 8,
        alpha = 0.8,
    )
    scatter!(ax, points_upstream; markersize, color = upstream_color, label = L"n_\text{u,sf}")

    scatter!(ax, points_downstream; markersize, color = downstream_color, label = L"n_\text{d,sf}")

    # velocity vectors
    arrows2d!(
        ax, v_downstream, fill(Point2(u_downstream, 0), length(v_downstream));
        color = downstream_color,
        # label = L"u_\text{d,sf}",
        shaftwidth,
        tiplength = 4, tipwidth = 8,
        alpha = 0.8,
    )

    vlines!(ax, 0, linestyle = :dashdot, color = :black)

    # Legend(fig[0,1], ax, merge=true, orientation = :horizontal, framevisible = false)

    fig
end

# ╔═╡ 8ad73454-040e-46da-8a83-cc3e467c2996
save("schematic.svg", schematic_figure)

# ╔═╡ 0a8a39f3-619d-402a-b6f2-4a4a4a11e794
save("schematic.png", schematic_figure, px_per_unit=2)

# ╔═╡ Cell order:
# ╠═67d360be-2973-11f1-85b4-13446ccb6184
# ╠═67d360be-2973-11f1-81b0-e9f876bff0c9
# ╠═e2eddf74-c76a-43e0-8c44-d1b06894c3cc
# ╠═2d19b46a-7de8-4156-b47b-1292642b5a54
# ╠═7d281913-c40d-4452-941f-7f5a28d62328
# ╠═df4b7c63-cb3e-4432-a4db-78f14f093a06
# ╠═6030ac76-e906-46fa-a68f-20a41c7ff130
# ╠═67d360be-2973-11f1-9cd7-4d87cb3d3c08
# ╠═8cf4fc55-daac-4915-ab30-6150fd948efd
# ╠═f2ac1dee-3b4a-4f75-b39e-5e29c457a9b2
# ╠═808a2a5f-a6e3-4660-9118-c5a430d20c34
# ╠═ce51675e-210d-4626-99fd-febf01108d67
# ╠═15823133-b986-4d03-be43-1d01a4aa74da
# ╠═041d4f92-c2e5-435a-a07f-31339554583e
# ╠═dd966d52-0d58-4ce2-a1d2-a3d2d0e45662
# ╠═1245e54f-c85c-445d-956c-117a01115126
# ╠═928a7e9c-9010-4eb4-87c9-10213f7a3164
# ╠═218cf991-e152-4b45-ae84-ee2d5d686e28
# ╠═1b47fa68-d9d1-486d-88dd-2c094248b974
# ╠═bbc09185-bb2b-4bf0-84fc-3d644851c3f3
# ╠═8ad73454-040e-46da-8a83-cc3e467c2996
# ╠═0a8a39f3-619d-402a-b6f2-4a4a4a11e794
# ╠═c6a99716-b7f9-4586-b0e8-257ded2564de
# ╠═36f34b4f-a3f3-4637-9db6-e6d435e54d18
