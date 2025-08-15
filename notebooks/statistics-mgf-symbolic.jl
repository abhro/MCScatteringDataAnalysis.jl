### A Pluto.jl notebook ###
# v0.20.15

using Markdown
using InteractiveUtils

# ╔═╡ a877d102-099f-11f0-37d1-39f0ef1f87e3
import Pkg; Pkg.activate(Base.current_project(), io=devnull)

# ╔═╡ a357df88-fae3-40ad-b636-a3fbfd5358bb
using PlutoUI

# ╔═╡ c875992e-a423-4376-9c58-88fc7a7c1a8c
using Symbolics

# ╔═╡ 4bcb0ee6-a5f7-457e-8457-2cbc07bb40db
using BiNormalDistributions

# ╔═╡ ecc7d76f-f828-49ea-b15b-69a8865be999
using Distributions

# ╔═╡ 857d8253-936d-4226-b9af-87553e65e0fb
md"""
# Preamble
"""

# ╔═╡ c7a48d7b-4bf7-429e-89de-18bbf5a8a7fd
md"""
## Import packages
"""

# ╔═╡ b6466bc4-ff63-4d88-9913-c2b1b80ea9eb
md"""
## Configure notebook appearance
"""

# ╔═╡ 1861682e-b53c-4acb-88b6-7a20eaa5a43c
TableOfContents(depth = 6)

# ╔═╡ f1c4bf46-c020-4bae-a805-0e834c90c82d
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

# ╔═╡ d5fcbe3c-ff12-4394-a225-16d211024693
md"""
## Set up base variables
"""

# ╔═╡ 2d52e729-1c62-49a3-8acc-a9e19a3c2797
@syms λ μ₁ μ₂ σ₁ σ₂ t μ σ

# ╔═╡ 68746676-7121-4ef4-afbf-ded292130da6
const Dt = Differential(t)

# ╔═╡ 30de8ac9-3935-4e0a-b287-94884d8d30e8
md"""
# Background
"""

# ╔═╡ 77dc30bd-0e8a-4b2c-8574-3bcba5947d65
md"""
The central moment function of order _k_ is
```math
\operatorname{E}{\!\left[\left(x - μ_x\right)^k\right]}
```
where ``μ_x = \operatorname{E}[x]`` is the mean of the distribution.
"""

# ╔═╡ 845e8631-258d-48bf-94c6-e3e84591f66b
md"""
# Normal distribution
"""

# ╔═╡ 8fbeb891-3c41-4c28-a5e6-c8053be05c90
md"""
## Closed form central moment function
"""

# ╔═╡ 296ed25b-738f-4867-a70b-45a98f3b609c
md"""
For a normal distribution ``\mathcal{N}(μ, σ)``, the central moment generating function is
```math
\operatorname{E}{\!\left[\left(x - μ\right)^k\right]} = \begin{cases}
    0 & \text{if $k$ is odd} \\
    σ^k \left(k-1\right)!! & \text{if $k$ is even}
\end{cases}
```
"""

# ╔═╡ 769de7d8-cd30-4670-a146-8aa981ea0666
doublefactorial(n) = prod(n:-2:1);

# ╔═╡ cfd61008-5a1b-473a-944d-793dd273877b
cm(::Type{Normal}, j::Integer) = isodd(j) ? 0 : doublefactorial(j - 1) * σ^j;

# ╔═╡ 684c95a2-f4a0-4571-b03f-e0057fb03459
md"""
## Moment generating function
"""

# ╔═╡ 5e31d1e0-f779-44cd-b3f1-79110ee7845f
md"""
Moment generating function of normal distribution
"""

# ╔═╡ c2239cf1-9068-4ca1-9497-06a474b315cc
mgf_N = exp(t * μ + 1//2 * t^2 * σ^2)

# ╔═╡ 7ade5063-2727-426d-bc6f-c0209259f3e8
md"""
## Central moment generating function
"""

# ╔═╡ 65ff6f13-f6fb-48a7-84b1-600a0a45eabc
md"""
Central moment generating function of normal distribution:
"""

# ╔═╡ 43240130-8613-47b8-8b27-ce249b718359
cmgf_N = exp(-μ * t) * mgf_N |> simplify

# ╔═╡ ca77433f-3259-442e-aafa-5508479618dc
md"""
## List of moments
"""

# ╔═╡ ace37383-6410-48a7-809d-4cd3bc8ca344
md"""
### Second moment
"""

# ╔═╡ 5ad34919-ca13-40d1-b367-7327dc721c24
Dt²mgf_N = (Dt^2)(mgf_N) |> expand_derivatives |> simplify;

# ╔═╡ bb5ce3b1-2eb5-43ab-b6c3-4d086f148afc
Dt²cmgf_N = (Dt^2)(cmgf_N) |> expand_derivatives |> simplify;

# ╔═╡ a1ce3165-4a28-49b1-996b-8fc222097294
md"""
Second moment of normal distribution around 0:
"""

# ╔═╡ 78b839d0-275b-4c29-9054-e4e68df8958a
substitute(Dt²mgf_N, t=>0)

# ╔═╡ 30ba0c5a-c448-44bd-8b4a-d442e3cc8fd6
md"""
Second moment of normal distribution around mean:
"""

# ╔═╡ 1d6f854e-4a02-4ba0-b73c-ef87a3e9382b
substitute(Dt²cmgf_N, t=>0)

# ╔═╡ bdcab3cd-27b3-4742-96a8-da539339c41c
md"""
### Third moment
"""

# ╔═╡ 233439cb-1829-4d8f-aa67-32b496483738
(Dt^3)(mgf_N) |> expand_derivatives |> simplify;

# ╔═╡ 72ff20d5-612a-42e2-8181-295a11c0fe1b
Dt³mgf_N = mgf_N * (μ + t * σ^2) * (3σ^2 + (μ + t * σ^2)^2);

# ╔═╡ 85dcd819-741f-4529-b38e-e3e5c4fd5eb8
Dt³cmgf_N = (Dt^3)(cmgf_N) |> expand_derivatives;

# ╔═╡ 3f0f7e44-2878-486c-9f35-f00049bdacaf
md"""
Third normal moment around 0:
"""

# ╔═╡ 2a4d89fc-0e70-4113-8a22-ff2da608f32b
substitute(Dt³mgf_N, t=>0)

# ╔═╡ 950b5f65-e3b6-47f3-9740-aedf2a02ee9f
md"""
Third normal moment around mean:
"""

# ╔═╡ 33b9f1ed-6c14-405d-b737-df2b3c7f13eb
substitute(Dt³cmgf_N, t=>0)

# ╔═╡ 7bb99412-4ec4-4945-9ca8-b7e3027f3068
md"""
### Fourth moment
"""

# ╔═╡ 16371045-d5de-4aea-a382-128b8ac19d8e
(Dt^4)(mgf_N) |> expand_derivatives;

# ╔═╡ f6c5e581-d42e-4143-b2ee-9b22b81d171f
Dt⁴mgf_N = mgf_N * (3σ^4 + 6σ^2 * (μ + t * σ^2)^2 + (μ + t * σ^2)^4);

# ╔═╡ 0c87fcf9-21aa-4897-b933-68126e4b02da
md"""
Fourth normal moment around 0:
"""

# ╔═╡ 00f1ba0d-a0d7-47c9-b063-ecec526383d5
substitute(Dt⁴mgf_N, t=>0)

# ╔═╡ 1163ee1e-1603-472d-9dea-dd0c6fb7d13a
md"""
# Binormal distribution
"""

# ╔═╡ c4bb8639-480c-4061-ace6-5ef57cc20973
md"""
## Closed form central moment function
"""

# ╔═╡ 5d6e5af0-6bd1-4456-a665-127afbaf3557
md"""
## Moment generating function
"""

# ╔═╡ e193653b-a53d-4978-9381-40919e2b1f5d
md"""
mgf of binormal distribution
"""

# ╔═╡ b71edbb2-a0be-45b9-b62f-6fe78113b2a1
md"""
auxiliary variables
"""

# ╔═╡ 0296cbe4-9df5-4d23-8666-f00663be923d
N₁_subst = Dict(μ=>μ₁, σ=>σ₁);

# ╔═╡ 071b8892-b132-40b6-a238-b1a796802c13
N₂_subst = Dict(μ=>μ₂, σ=>σ₂);

# ╔═╡ b7c6e76d-55c1-4e33-a540-dcb57c6f4436
function cm(::Type{BiNormal}, j::Integer)
    moment = zero(λ)
    for k in 0:j
        normal_central = cm(Normal, k)
        moment += (
                 λ  * binomial(j, k) * (μ₁ - μ)^(j-k) * substitute(normal_central, N₁_subst)
            + (1-λ) * binomial(j, k) * (μ₂ - μ)^(j-k) * substitute(normal_central, N₂_subst)
        )
    end
    return moment
end

# ╔═╡ ad80264e-7851-4da6-9aa5-5739011fcbb6
cm.(Normal, 0:2:14)

# ╔═╡ 6deb3f39-cbb9-4330-a602-27fd88e89a1c
cm.(Normal, 1:2:15)

# ╔═╡ 58ae6e09-fdcc-408c-84d4-6f146b7a6a2c
mgf_BN = λ * substitute(mgf_N, N₁_subst) + (1 - λ) * substitute(mgf_N, N₂_subst)

# ╔═╡ 20dd0bc7-9e56-4db3-917a-8700a0fd9eef
md"""
## First moment
"""

# ╔═╡ 6547815e-53f8-4092-b50c-4eba80d5366c
DtM_BN = Dt(mgf_BN) |> expand_derivatives;

# ╔═╡ 19a1c7da-b60e-4d16-9c59-10e98540fbd6
const mean_binormal = substitute(DtM_BN, t=>0)

# ╔═╡ 5b20f228-e334-44f3-8ac2-c5c6d8c6dc20
md"""
## Central moment generating function
"""

# ╔═╡ b32e6a33-af1a-4991-9d54-339dad07b3f9
md"""
Central moment generating function of the bi-normal distribution:
"""

# ╔═╡ 0535663e-5723-4a46-b63d-8e371553cdea
cmgf_BN = exp(-mean_binormal * t) * mgf_BN

# ╔═╡ d7a30a6e-f523-43a2-9c05-7f6f2b147d8d
md"""
## List of higher-order moments
"""

# ╔═╡ 12bece90-4ead-41b1-aaf6-0e4f9a5cc81a
md"""
### Second moment
"""

# ╔═╡ 4572a3e6-71e4-4e0c-945d-127f2c93427f
Dt²mgf_BN = (Dt^2)(mgf_BN) |> expand_derivatives;

# ╔═╡ f07c41e4-f694-4259-ae1b-1d9fa407cd74
Dt²cmgf_BN = (Dt^2)(cmgf_BN) |> expand_derivatives;

# ╔═╡ 043a1644-437c-479e-afaa-791c3d8f81b7
md"""
Second moment of binormal around 0:
"""

# ╔═╡ a685362b-a0ff-4f08-9a1f-6b05d34fb23d
substitute(Dt²mgf_BN, t => 0) |> simplify

# ╔═╡ 602f7abd-ad30-4518-a9fb-9c1d5a1c6867
md"""
Second moment of binormal around mean:
"""

# ╔═╡ 9448227f-8c5d-4535-80a6-eb5bd4c31c59
substitute(Dt²cmgf_BN, t => 0) |> simplify

# ╔═╡ a646c9ca-66d8-4dc9-bc16-49cec6fe115d
substitute(cm(BiNormal, 2), μ => mean_binormal) |> simplify

# ╔═╡ 600915b5-eb79-456d-b7b8-debc878b6c74
variance_BN = λ*σ₁^2 + (1-λ)*σ₂^2 + λ*(1-λ)*(μ₁-μ₂)^2

# ╔═╡ c3674e3c-f623-4590-8b62-349304055688
md"""
### Third moment
"""

# ╔═╡ e24afe40-4672-49bc-968c-b042a2fa604a
Dt³mgf_BN = λ * substitute(Dt³mgf_N, N₁_subst) + (1-λ) * substitute(Dt³mgf_N, N₂_subst);

# ╔═╡ 94726946-6761-4e63-a59e-9b5b36082eed
Dt³cmgf_BN = (Dt^3)(cmgf_BN) |> expand_derivatives |> simplify;

# ╔═╡ 923a3ec9-614f-4283-8f58-4216fa175d23
md"""
Third binormal moment around 0
"""

# ╔═╡ a9748d7b-a604-4faa-934a-b10e43e04966
substitute(Dt³mgf_BN, t => 0) |> simplify

# ╔═╡ c2adfec5-00a0-405d-aa02-5999f727aba6
md"""
Third binormal moment around mean
"""

# ╔═╡ 0fd3e6e4-7cf7-4c7a-8df2-6de722d1476f
(substitute(Dt³cmgf_BN, t => 0) |> simplify |> terms)[[2,3,7,1,4,5,6]]

# ╔═╡ a396f0b6-4630-4ef2-b428-4fa8f2d8a593
cm(BiNormal, 3)

# ╔═╡ e7c35714-a530-4d39-82a1-4bac1070e7ff
cm3_BN = (
         λ  * μ₁^3 + 3 *    λ  * μ₁ * σ₁^2
    + (1-λ) * μ₂^3 + 3 * (1-λ) * μ₂ * σ₁^2
    + 2 * (λ*μ₁ + (1-λ)*μ₂)^3
    - 3 * (λ*μ₁ + (1-λ)*μ₂) * (λ*(μ₁^2 + σ₁^2) + (1-λ)*(μ₂^2 + σ₂^2))
)

# ╔═╡ 56d4d179-7c5d-4918-a8ab-41a330a0c183
md"""
#### Skewness
"""

# ╔═╡ 97ec97c1-b812-4e28-809d-ece86f94bd59
md"""
skewness = (third central moment) / (second central moment)^(3/2)
"""

# ╔═╡ bd9b02ae-777c-4147-aecd-065de1889471
md"""
see the docs for BiNormalDistributions.jl
"""

# ╔═╡ c04ee9dc-dc26-469e-801b-ea8f1cc251a2
md"""
### Fourth moment
"""

# ╔═╡ f239f933-d8bc-44de-912a-fedf12101625
((Dt^4)(mgf_BN) |> expand_derivatives |> simplify |> terms);#[[1,2,4,6,5,3]]

# ╔═╡ a4625d7e-6051-4fd8-9f85-866816760d04
Dt⁴mgf_BN = (
    λ * exp(t*μ₁ + 1//2 * t^2*σ₁^2) * (
        (μ₁ + t * σ₁^2)^4 + 3σ₁^4 + 6σ₁^2 * (μ₁ + t*σ₁^2)^2)
    +
    (1-λ) * exp(t*μ₂ + 1//2 * t^2*σ₂^2) * (
        (μ₂ + t * σ₂^2)^4 + 3σ₂^4 + 6σ₂^2 * (μ₂ + t*σ₂^2)^2)
);

# ╔═╡ 72a2e24a-3285-451f-b68a-40936067a84f
Dt⁴cmgf_BN = (Dt^4)(cmgf_BN) |> expand_derivatives |> simplify;

# ╔═╡ 6067ded8-26e5-4a3a-99d3-e13c1ed66484
md"""
Fourth binormal moment around 0:
"""

# ╔═╡ 91a9a42d-b323-4916-a025-87aba0988d4e
substitute(Dt⁴mgf_BN, t=>0)

# ╔═╡ 09e9c2e5-c2d1-43cd-b9f0-8474dd4d3ca7
md"""
Fourth binormal moment around mean:
"""

# ╔═╡ d6364bcd-dee0-481c-a61a-48c65654b0c1
cm4_BN = substitute(Dt⁴cmgf_BN, t=>0) |> simplify

# ╔═╡ cce3d179-b5c9-4124-ac83-66348568cd60
cm4_BN_terms = terms(cm4_BN)

# ╔═╡ 30ddf59e-8714-4c78-9acb-4d5fbd1f54e9
cm4_BN_terms[[7,1,5,10,3,2,9,4,8,6]]

# ╔═╡ a5d62918-69a2-48ee-ad5d-760d9546f475
substitute(cm(BiNormal, 4), μ=>mean_binormal) |> simplify |> terms

# ╔═╡ c210894e-e087-4f89-ac84-516f89e919c6
md"""
#### Kurtosis
"""

# ╔═╡ 57a7711d-99a1-4b2f-9b24-d1b26ed1131d
md"""
kurtosis = (fourth central moment) / (second central moment)^2
"""

# ╔═╡ f49ae510-2b40-4b16-aee2-b6815582abc2
md"""
### Fifth moment
"""

# ╔═╡ 73195df0-ad44-43e9-8366-fc66707f3a30
Dt⁵mgf_BN = (Dt^5)(mgf_BN) |> expand_derivatives |> simplify

# ╔═╡ e8341209-c5b1-4298-88e6-7f7b19d94dc6
Dt⁵cmgf_BN = (Dt^5)(cmgf_BN) |> expand_derivatives |> simplify

# ╔═╡ 7b774191-d658-4b64-8a7a-ea9b940feb17
substitute(Dt⁵mgf_BN, t=> 0)

# ╔═╡ 32170b91-c863-4817-8c55-5fbe37024872
cm5_BN = substitute(Dt⁵cmgf_BN, t=> 0) |> simplify

# ╔═╡ 5f833a90-3b4b-43e3-9c43-1dfae492f3db
cm5_BN |> terms

# ╔═╡ 171bc4f4-0760-44ac-bf99-066e0a26eef9
substitute(cm(BiNormal, 5), μ=>mean_binormal) |> simplify |> terms

# ╔═╡ 5e92bf1d-5481-495b-96e0-2a3096f05468
md"""
### Sixth moment
"""

# ╔═╡ 39d814cf-d038-40e7-af90-4066c987eb21
md"""
Sixth moment around mean
"""

# ╔═╡ de52b701-ead7-48c0-a60f-c58384e846b9
cm(BiNormal, 6)

# ╔═╡ 9cadf5aa-bb63-4179-8a99-04714eb74617
substitute(cm(BiNormal, 6), μ=>mean_binormal) |> simplify |> terms

# ╔═╡ Cell order:
# ╟─857d8253-936d-4226-b9af-87553e65e0fb
# ╟─c7a48d7b-4bf7-429e-89de-18bbf5a8a7fd
# ╠═a877d102-099f-11f0-37d1-39f0ef1f87e3
# ╠═a357df88-fae3-40ad-b636-a3fbfd5358bb
# ╠═c875992e-a423-4376-9c58-88fc7a7c1a8c
# ╟─b6466bc4-ff63-4d88-9913-c2b1b80ea9eb
# ╠═1861682e-b53c-4acb-88b6-7a20eaa5a43c
# ╟─f1c4bf46-c020-4bae-a805-0e834c90c82d
# ╟─d5fcbe3c-ff12-4394-a225-16d211024693
# ╠═2d52e729-1c62-49a3-8acc-a9e19a3c2797
# ╠═68746676-7121-4ef4-afbf-ded292130da6
# ╟─30de8ac9-3935-4e0a-b287-94884d8d30e8
# ╟─77dc30bd-0e8a-4b2c-8574-3bcba5947d65
# ╟─845e8631-258d-48bf-94c6-e3e84591f66b
# ╟─8fbeb891-3c41-4c28-a5e6-c8053be05c90
# ╟─296ed25b-738f-4867-a70b-45a98f3b609c
# ╠═769de7d8-cd30-4670-a146-8aa981ea0666
# ╠═cfd61008-5a1b-473a-944d-793dd273877b
# ╠═ad80264e-7851-4da6-9aa5-5739011fcbb6
# ╠═6deb3f39-cbb9-4330-a602-27fd88e89a1c
# ╟─684c95a2-f4a0-4571-b03f-e0057fb03459
# ╟─5e31d1e0-f779-44cd-b3f1-79110ee7845f
# ╟─c2239cf1-9068-4ca1-9497-06a474b315cc
# ╟─7ade5063-2727-426d-bc6f-c0209259f3e8
# ╟─65ff6f13-f6fb-48a7-84b1-600a0a45eabc
# ╟─43240130-8613-47b8-8b27-ce249b718359
# ╟─ca77433f-3259-442e-aafa-5508479618dc
# ╟─ace37383-6410-48a7-809d-4cd3bc8ca344
# ╠═5ad34919-ca13-40d1-b367-7327dc721c24
# ╠═bb5ce3b1-2eb5-43ab-b6c3-4d086f148afc
# ╟─a1ce3165-4a28-49b1-996b-8fc222097294
# ╟─78b839d0-275b-4c29-9054-e4e68df8958a
# ╟─30ba0c5a-c448-44bd-8b4a-d442e3cc8fd6
# ╟─1d6f854e-4a02-4ba0-b73c-ef87a3e9382b
# ╟─bdcab3cd-27b3-4742-96a8-da539339c41c
# ╠═233439cb-1829-4d8f-aa67-32b496483738
# ╠═72ff20d5-612a-42e2-8181-295a11c0fe1b
# ╠═85dcd819-741f-4529-b38e-e3e5c4fd5eb8
# ╟─3f0f7e44-2878-486c-9f35-f00049bdacaf
# ╟─2a4d89fc-0e70-4113-8a22-ff2da608f32b
# ╟─950b5f65-e3b6-47f3-9740-aedf2a02ee9f
# ╟─33b9f1ed-6c14-405d-b737-df2b3c7f13eb
# ╟─7bb99412-4ec4-4945-9ca8-b7e3027f3068
# ╠═16371045-d5de-4aea-a382-128b8ac19d8e
# ╠═f6c5e581-d42e-4143-b2ee-9b22b81d171f
# ╟─0c87fcf9-21aa-4897-b933-68126e4b02da
# ╠═00f1ba0d-a0d7-47c9-b063-ecec526383d5
# ╟─1163ee1e-1603-472d-9dea-dd0c6fb7d13a
# ╠═4bcb0ee6-a5f7-457e-8457-2cbc07bb40db
# ╠═ecc7d76f-f828-49ea-b15b-69a8865be999
# ╟─c4bb8639-480c-4061-ace6-5ef57cc20973
# ╠═b7c6e76d-55c1-4e33-a540-dcb57c6f4436
# ╟─5d6e5af0-6bd1-4456-a665-127afbaf3557
# ╟─e193653b-a53d-4978-9381-40919e2b1f5d
# ╟─58ae6e09-fdcc-408c-84d4-6f146b7a6a2c
# ╟─b71edbb2-a0be-45b9-b62f-6fe78113b2a1
# ╠═0296cbe4-9df5-4d23-8666-f00663be923d
# ╠═071b8892-b132-40b6-a238-b1a796802c13
# ╟─20dd0bc7-9e56-4db3-917a-8700a0fd9eef
# ╠═6547815e-53f8-4092-b50c-4eba80d5366c
# ╠═19a1c7da-b60e-4d16-9c59-10e98540fbd6
# ╟─5b20f228-e334-44f3-8ac2-c5c6d8c6dc20
# ╟─b32e6a33-af1a-4991-9d54-339dad07b3f9
# ╠═0535663e-5723-4a46-b63d-8e371553cdea
# ╟─d7a30a6e-f523-43a2-9c05-7f6f2b147d8d
# ╟─12bece90-4ead-41b1-aaf6-0e4f9a5cc81a
# ╠═4572a3e6-71e4-4e0c-945d-127f2c93427f
# ╠═f07c41e4-f694-4259-ae1b-1d9fa407cd74
# ╟─043a1644-437c-479e-afaa-791c3d8f81b7
# ╠═a685362b-a0ff-4f08-9a1f-6b05d34fb23d
# ╟─602f7abd-ad30-4518-a9fb-9c1d5a1c6867
# ╠═9448227f-8c5d-4535-80a6-eb5bd4c31c59
# ╠═a646c9ca-66d8-4dc9-bc16-49cec6fe115d
# ╠═600915b5-eb79-456d-b7b8-debc878b6c74
# ╟─c3674e3c-f623-4590-8b62-349304055688
# ╠═e24afe40-4672-49bc-968c-b042a2fa604a
# ╠═94726946-6761-4e63-a59e-9b5b36082eed
# ╟─923a3ec9-614f-4283-8f58-4216fa175d23
# ╠═a9748d7b-a604-4faa-934a-b10e43e04966
# ╟─c2adfec5-00a0-405d-aa02-5999f727aba6
# ╠═0fd3e6e4-7cf7-4c7a-8df2-6de722d1476f
# ╠═a396f0b6-4630-4ef2-b428-4fa8f2d8a593
# ╠═e7c35714-a530-4d39-82a1-4bac1070e7ff
# ╟─56d4d179-7c5d-4918-a8ab-41a330a0c183
# ╟─97ec97c1-b812-4e28-809d-ece86f94bd59
# ╟─bd9b02ae-777c-4147-aecd-065de1889471
# ╟─c04ee9dc-dc26-469e-801b-ea8f1cc251a2
# ╠═f239f933-d8bc-44de-912a-fedf12101625
# ╠═a4625d7e-6051-4fd8-9f85-866816760d04
# ╠═72a2e24a-3285-451f-b68a-40936067a84f
# ╟─6067ded8-26e5-4a3a-99d3-e13c1ed66484
# ╠═91a9a42d-b323-4916-a025-87aba0988d4e
# ╟─09e9c2e5-c2d1-43cd-b9f0-8474dd4d3ca7
# ╠═d6364bcd-dee0-481c-a61a-48c65654b0c1
# ╠═cce3d179-b5c9-4124-ac83-66348568cd60
# ╠═30ddf59e-8714-4c78-9acb-4d5fbd1f54e9
# ╠═a5d62918-69a2-48ee-ad5d-760d9546f475
# ╟─c210894e-e087-4f89-ac84-516f89e919c6
# ╟─57a7711d-99a1-4b2f-9b24-d1b26ed1131d
# ╟─f49ae510-2b40-4b16-aee2-b6815582abc2
# ╠═73195df0-ad44-43e9-8366-fc66707f3a30
# ╠═e8341209-c5b1-4298-88e6-7f7b19d94dc6
# ╠═7b774191-d658-4b64-8a7a-ea9b940feb17
# ╠═32170b91-c863-4817-8c55-5fbe37024872
# ╠═5f833a90-3b4b-43e3-9c43-1dfae492f3db
# ╠═171bc4f4-0760-44ac-bf99-066e0a26eef9
# ╟─5e92bf1d-5481-495b-96e0-2a3096f05468
# ╟─39d814cf-d038-40e7-af90-4066c987eb21
# ╠═de52b701-ead7-48c0-a60f-c58384e846b9
# ╠═9cadf5aa-bb63-4179-8a99-04714eb74617
