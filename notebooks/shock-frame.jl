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
using CairoMakie; set_theme!(theme_latexfonts())

# ╔═╡ e2eddf74-c76a-43e0-8c44-d1b06894c3cc
using Distributions

# ╔═╡ 2d19b46a-7de8-4156-b47b-1292642b5a54
using Random

# ╔═╡ 9b4eccbe-2f83-4fd3-99fa-cccb192fd270
using Sobol

# ╔═╡ 7d281913-c40d-4452-941f-7f5a28d62328
rng = Xoshiro(5)

# ╔═╡ df4b7c63-cb3e-4432-a4db-78f14f093a06
# compression ratio
r = 7;

# ╔═╡ 6030ac76-e906-46fa-a68f-20a41c7ff130
# number of upstream particles
n_upstream = 50;

# ╔═╡ 67d360be-2973-11f1-9cd7-4d87cb3d3c08
# number of downstream_particles
n_downstream = n_upstream * r

# ╔═╡ 8cf4fc55-daac-4915-ab30-6150fd948efd
y_lo, y_hi = -2, 2;

# ╔═╡ d6818676-dd60-4fe4-813b-268a849e9cb5
points_downstream = let
    x_downstream_lo = 0.1
    x_downstream_hi = 4.4
    points_downstream_seq = SobolSeq([x_downstream_lo, y_lo], [x_downstream_hi, y_hi])
    map(Point2, next!(points_downstream_seq) for i in 1:n_downstream)
end

# ╔═╡ 15823133-b986-4d03-be43-1d01a4aa74da
points_upstream = let
    x_upstream_lo = -5
    x_upstream_hi = -0.1
    points_upstream_seq = SobolSeq([x_upstream_lo, y_lo], [x_upstream_hi, y_hi])
    map(Point2, next!(points_upstream_seq) for i in 1:n_upstream)
end

# ╔═╡ 655afaea-c38d-49b3-9897-f98168e14036
# number of upstream particles that get a velocity vector
# vn_upstream = 40
vn_upstream = n_upstream ÷ 2

# ╔═╡ 83065981-47b3-42ea-9f61-079abe92e6f6
# number of downstream particles that get a velocity vector
# vn_downstream = 90
vn_downstream = n_downstream ÷ (8//5)

# ╔═╡ 041d4f92-c2e5-435a-a07f-31339554583e
begin
    # downstream particles that get a velocity vector
    v_downstream_points = sample(rng, points_downstream, vn_downstream, replace = false)
    # upstream particles that get a velocity vector
    v_upstream_points = sample(rng, points_upstream, vn_upstream, replace = false)
end;

# ╔═╡ dd966d52-0d58-4ce2-a1d2-a3d2d0e45662
# size of arrow for upstream velocity
u_upstream = 0.45 * 2

# ╔═╡ 1245e54f-c85c-445d-956c-117a01115126
# size of arrow for downstream velocity
u_downstream = u_upstream / r^0.8

# ╔═╡ 8e6765b5-04c6-4111-aee4-6af723954014
#upstream_velocities = fill(Point2(u_upstream, 0), length(v_upstream))
upstream_velocities = let
    # small fluctuations in angle due to low upstream temperature
    fluctuation_dist = Uniform(-1/150, 1/150) # in units of π (straight angle)
    fluctuation_angles = rand(rng, fluctuation_dist, vn_upstream)
    # Makie has the y-axis as 0 degrees, so add a right-angle clockwise
    fluctuations = Point2.(sincospi.(1/2 .+ fluctuation_angles))
    u_upstream .* fluctuations
end

# ╔═╡ 01e587a5-50a7-4990-bb3e-32f6d301ff0e
downstream_velocities_thermal = let
    # at high downstream temperature, the entire unit circle is fair game
    θ_dist = Uniform(0, 2π)
    fluctuation_angles = rand(rng, θ_dist, vn_downstream)
    # unlike upstream, all angles are likely, so no right-angle offset needed
    fluctuations = Point2.(sincos.(fluctuation_angles))
    u_downstream/1.3 .* fluctuations
end

# ╔═╡ 5ae39433-ee6b-413e-9217-2fa200e9c012
downstream_velocities_bulk = fill(Point2(u_downstream, 0), vn_downstream)

# ╔═╡ e87b03ea-0fca-49d2-a3aa-6185d9c5c721
downstream_velocities = downstream_velocities_bulk + downstream_velocities_thermal

# ╔═╡ 928a7e9c-9010-4eb4-87c9-10213f7a3164
markersize = 5;

# ╔═╡ 218cf991-e152-4b45-ae84-ee2d5d686e28
# plotting properties for velocity arrows
arrow_properties = (; shaftwidth = 1, tiplength = 4, tipwidth = 8, alpha = 0.6);

# ╔═╡ 9388a9d7-38b7-4cf2-b1c4-52a10aef70ef
upstream_text = L"\textbf{upstream}, \; u_\text{u,sf}, \; n_\text{u,sf}, \; T_\text{u,sf}"

# ╔═╡ 41af110d-505f-451c-9992-36f20f8e4e01
downstream_text = L"\textbf{downstream}, \; u_\text{d,sf}, \; n_\text{d,sf}, \; T_\text{d,sf}"

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

    text!(ax, -0.7, 2.2, text = L"\textbf{shock front}")
    text!(ax, -4.6, 2.2, text = upstream_text,  color = upstream_color)
    text!(ax, 1.4, 2.2, text = downstream_text, color = downstream_color)

    # plot the particle populations as dots
    scatter!(ax, points_upstream; markersize, color = upstream_color, label = L"n_\text{u,sf}")
    scatter!(ax, points_downstream; markersize, color = downstream_color, label = L"n_\text{d,sf}")

    # attach velocity vector arrows to the particles
    arrows2d!(
        ax, v_upstream_points, upstream_velocities;
        color = upstream_color,
        # label = L"u_\text{u,sf}",
        arrow_properties...,
    )
    arrows2d!(
        ax, v_downstream_points, downstream_velocities;
        color = downstream_color,
        # label = L"u_\text{d,sf}",
        arrow_properties...,
    )

    vlines!(ax, 0, ymax = 0.95, linestyle = :dashdot, color = :black)

    # Legend(fig[0,1], ax, merge=true, orientation = :horizontal, framevisible = false)

    fig
end

# ╔═╡ 8ad73454-040e-46da-8a83-cc3e467c2996
# ╠═╡ disabled = true
#=╠═╡
begin
    save("shock-frame-flow-schematic.svg", schematic_figure)
    save("shock-frame-flow-schematic.png", schematic_figure, px_per_unit=2)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═67d360be-2973-11f1-85b4-13446ccb6184
# ╠═67d360be-2973-11f1-81b0-e9f876bff0c9
# ╠═e2eddf74-c76a-43e0-8c44-d1b06894c3cc
# ╠═2d19b46a-7de8-4156-b47b-1292642b5a54
# ╠═9b4eccbe-2f83-4fd3-99fa-cccb192fd270
# ╠═7d281913-c40d-4452-941f-7f5a28d62328
# ╠═df4b7c63-cb3e-4432-a4db-78f14f093a06
# ╠═6030ac76-e906-46fa-a68f-20a41c7ff130
# ╠═67d360be-2973-11f1-9cd7-4d87cb3d3c08
# ╠═8cf4fc55-daac-4915-ab30-6150fd948efd
# ╠═d6818676-dd60-4fe4-813b-268a849e9cb5
# ╠═15823133-b986-4d03-be43-1d01a4aa74da
# ╠═655afaea-c38d-49b3-9897-f98168e14036
# ╠═83065981-47b3-42ea-9f61-079abe92e6f6
# ╠═041d4f92-c2e5-435a-a07f-31339554583e
# ╠═dd966d52-0d58-4ce2-a1d2-a3d2d0e45662
# ╠═1245e54f-c85c-445d-956c-117a01115126
# ╠═8e6765b5-04c6-4111-aee4-6af723954014
# ╠═01e587a5-50a7-4990-bb3e-32f6d301ff0e
# ╠═5ae39433-ee6b-413e-9217-2fa200e9c012
# ╠═e87b03ea-0fca-49d2-a3aa-6185d9c5c721
# ╠═928a7e9c-9010-4eb4-87c9-10213f7a3164
# ╠═218cf991-e152-4b45-ae84-ee2d5d686e28
# ╟─bbc09185-bb2b-4bf0-84fc-3d644851c3f3
# ╠═8ad73454-040e-46da-8a83-cc3e467c2996
# ╠═9388a9d7-38b7-4cf2-b1c4-52a10aef70ef
# ╠═41af110d-505f-451c-9992-36f20f8e4e01
# ╠═c6a99716-b7f9-4586-b0e8-257ded2564de
# ╠═36f34b4f-a3f3-4637-9db6-e6d435e54d18
