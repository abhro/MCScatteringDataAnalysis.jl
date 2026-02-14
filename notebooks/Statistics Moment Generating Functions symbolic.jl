### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ a357df88-fae3-40ad-b636-a3fbfd5358bb
using PlutoUI

# ╔═╡ c875992e-a423-4376-9c58-88fc7a7c1a8c
using Symbolics, Latexify

# ╔═╡ ecc7d76f-f828-49ea-b15b-69a8865be999
using Distributions

# ╔═╡ 589dae68-0c2d-4c5a-a58d-2d09f0f51c74
md"""
# Derivation of moments for the bi-normal distribution
"""

# ╔═╡ 857d8253-936d-4226-b9af-87553e65e0fb
md"""
## Preamble
"""

# ╔═╡ c7a48d7b-4bf7-429e-89de-18bbf5a8a7fd
md"""
### Import packages
"""

# ╔═╡ e8d2ef67-e395-4915-9b19-4df3e10327dc
md"""
The following two cells relate to the use of DrWatson.jl. They've been deactivated since the project's dependency list in Project.toml no longer contains the packages used in this notebook for performance reasons. If they're added back, feel free to re-enable the cells.
"""

# ╔═╡ a877d102-099f-11f0-37d1-39f0ef1f87e3
# using DrWatson

# ╔═╡ e5e0e4e2-2df1-4536-9cc5-bdcec6fc13de
# @quickactivate "MCScatteringDataAnalysis"

# ╔═╡ b6466bc4-ff63-4d88-9913-c2b1b80ea9eb
md"""
### Configure notebook appearance
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
### Set up base variables
"""

# ╔═╡ 2d52e729-1c62-49a3-8acc-a9e19a3c2797
@syms λ μ₁ μ₂ σ₁ σ₂ t μ σ

# ╔═╡ 68746676-7121-4ef4-afbf-ded292130da6
const Dt = Differential(t)

# ╔═╡ 30de8ac9-3935-4e0a-b287-94884d8d30e8
md"""
## Background
"""

# ╔═╡ afc8be14-3e99-4eab-bbbd-673ec87e1052
md"""
The central moment function of order _k_ is
```math
\operatorname{E}{\!\left[x^k\right]}.
```
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
## Normal distribution
"""

# ╔═╡ 8fbeb891-3c41-4c28-a5e6-c8053be05c90
md"""
### Closed form central moment function
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
### Moment generating function
"""

# ╔═╡ 5e31d1e0-f779-44cd-b3f1-79110ee7845f
md"""
Moment generating function of normal distribution
"""

# ╔═╡ c2239cf1-9068-4ca1-9497-06a474b315cc
mgf_N = exp(t * μ + 1//2 * t^2 * σ^2)

# ╔═╡ 7ade5063-2727-426d-bc6f-c0209259f3e8
md"""
### Central moment generating function
"""

# ╔═╡ 65ff6f13-f6fb-48a7-84b1-600a0a45eabc
md"""
Central moment generating function of normal distribution:
"""

# ╔═╡ 43240130-8613-47b8-8b27-ce249b718359
cmgf_N = exp(-μ * t) * mgf_N |> simplify

# ╔═╡ ca77433f-3259-442e-aafa-5508479618dc
md"""
### List of moments
"""

# ╔═╡ ace37383-6410-48a7-809d-4cd3bc8ca344
md"""
#### Second moment
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
substitute(Dt²mgf_N, t => 0)

# ╔═╡ 30ba0c5a-c448-44bd-8b4a-d442e3cc8fd6
md"""
Second moment of normal distribution around mean:
"""

# ╔═╡ 1d6f854e-4a02-4ba0-b73c-ef87a3e9382b
substitute(Dt²cmgf_N, t => 0)

# ╔═╡ bdcab3cd-27b3-4742-96a8-da539339c41c
md"""
#### Third moment
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
substitute(Dt³mgf_N, t => 0)

# ╔═╡ 950b5f65-e3b6-47f3-9740-aedf2a02ee9f
md"""
Third normal moment around mean:
"""

# ╔═╡ 33b9f1ed-6c14-405d-b737-df2b3c7f13eb
substitute(Dt³cmgf_N, t => 0)

# ╔═╡ 7bb99412-4ec4-4945-9ca8-b7e3027f3068
md"""
#### Fourth moment
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
substitute(Dt⁴mgf_N, t => 0)

# ╔═╡ 1163ee1e-1603-472d-9dea-dd0c6fb7d13a
md"""
## Binormal distribution
"""

# ╔═╡ 4bcb0ee6-a5f7-457e-8457-2cbc07bb40db
# using BiNormalDistributions

# ╔═╡ c4bb8639-480c-4061-ace6-5ef57cc20973
md"""
### Closed form central moment function
"""

# ╔═╡ 5d6e5af0-6bd1-4456-a665-127afbaf3557
md"""
### Moment generating function
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
N₁_subst = Dict(μ => μ₁, σ => σ₁);

# ╔═╡ 071b8892-b132-40b6-a238-b1a796802c13
N₂_subst = Dict(μ => μ₂, σ => σ₂);

# ╔═╡ b7c6e76d-55c1-4e33-a540-dcb57c6f4436
# function cm(::Type{BiNormal}, j::Integer)
function cm(::Val{:BiNormal}, j::Integer)
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
### First moment
"""

# ╔═╡ 6547815e-53f8-4092-b50c-4eba80d5366c
DtM_BN = Dt(mgf_BN) |> expand_derivatives;

# ╔═╡ 19a1c7da-b60e-4d16-9c59-10e98540fbd6
const mean_binormal = substitute(DtM_BN, t => 0)

# ╔═╡ 5b20f228-e334-44f3-8ac2-c5c6d8c6dc20
md"""
### Central moment generating function
"""

# ╔═╡ b32e6a33-af1a-4991-9d54-339dad07b3f9
md"""
Central moment generating function of the bi-normal distribution:
"""

# ╔═╡ 0535663e-5723-4a46-b63d-8e371553cdea
cmgf_BN = exp(-mean_binormal * t) * mgf_BN

# ╔═╡ d7a30a6e-f523-43a2-9c05-7f6f2b147d8d
md"""
### List of higher-order moments
"""

# ╔═╡ 12bece90-4ead-41b1-aaf6-0e4f9a5cc81a
md"""
#### Second moment
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
# substitute(cm(BiNormal, 2), μ => mean_binormal) |> simplify
substitute(cm(Val(:BiNormal), 2), μ => mean_binormal) |> simplify

# ╔═╡ 600915b5-eb79-456d-b7b8-debc878b6c74
variance_BN = λ*σ₁^2 + (1-λ)*σ₂^2 + λ*(1-λ)*(μ₁-μ₂)^2

# ╔═╡ c3674e3c-f623-4590-8b62-349304055688
md"""
#### Third moment
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
# cm(BiNormal, 3)
cm(Val(:BiNormal), 3)

# ╔═╡ e7c35714-a530-4d39-82a1-4bac1070e7ff
cm3_BN = (
         λ  * μ₁^3 + 3 *    λ  * μ₁ * σ₁^2
    + (1-λ) * μ₂^3 + 3 * (1-λ) * μ₂ * σ₁^2
    + 2 * (λ*μ₁ + (1-λ)*μ₂)^3
    - 3 * (λ*μ₁ + (1-λ)*μ₂) * (λ*(μ₁^2 + σ₁^2) + (1-λ)*(μ₂^2 + σ₂^2))
)

# ╔═╡ 56d4d179-7c5d-4918-a8ab-41a330a0c183
md"""
##### Skewness
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
#### Fourth moment
"""

# ╔═╡ f239f933-d8bc-44de-912a-fedf12101625
((Dt^4)(mgf_BN) |> expand_derivatives |> simplify |> terms); #[[1,2,4,6,5,3]]

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
substitute(Dt⁴mgf_BN, t => 0)

# ╔═╡ 09e9c2e5-c2d1-43cd-b9f0-8474dd4d3ca7
md"""
Fourth binormal moment around mean:
"""

# ╔═╡ d6364bcd-dee0-481c-a61a-48c65654b0c1
cm4_BN = substitute(Dt⁴cmgf_BN, t => 0) |> simplify

# ╔═╡ cce3d179-b5c9-4124-ac83-66348568cd60
cm4_BN_terms = terms(cm4_BN)

# ╔═╡ a5d62918-69a2-48ee-ad5d-760d9546f475
substitute(cm(Val(:BiNormal), 4), μ => mean_binormal) |> simplify |> terms

# ╔═╡ c210894e-e087-4f89-ac84-516f89e919c6
md"""
##### Kurtosis
"""

# ╔═╡ 57a7711d-99a1-4b2f-9b24-d1b26ed1131d
md"""
kurtosis = (fourth central moment) / (second central moment)^2
"""

# ╔═╡ f49ae510-2b40-4b16-aee2-b6815582abc2
md"""
#### Fifth moment
"""

# ╔═╡ 73195df0-ad44-43e9-8366-fc66707f3a30
Dt⁵mgf_BN = (Dt^5)(mgf_BN) |> expand_derivatives |> simplify

# ╔═╡ e8341209-c5b1-4298-88e6-7f7b19d94dc6
Dt⁵cmgf_BN = (Dt^5)(cmgf_BN) |> expand_derivatives |> simplify

# ╔═╡ 7b774191-d658-4b64-8a7a-ea9b940feb17
substitute(Dt⁵mgf_BN, t => 0)

# ╔═╡ 32170b91-c863-4817-8c55-5fbe37024872
cm5_BN = substitute(Dt⁵cmgf_BN, t => 0) |> simplify

# ╔═╡ 5f833a90-3b4b-43e3-9c43-1dfae492f3db
cm5_BN |> terms

# ╔═╡ 171bc4f4-0760-44ac-bf99-066e0a26eef9
substitute(cm(Val(:BiNormal), 5), μ => mean_binormal) |> simplify |> terms

# ╔═╡ 5e92bf1d-5481-495b-96e0-2a3096f05468
md"""
#### Sixth moment
"""

# ╔═╡ 39d814cf-d038-40e7-af90-4066c987eb21
md"""
Sixth moment around mean
"""

# ╔═╡ de52b701-ead7-48c0-a60f-c58384e846b9
cm(Val(:BiNormal), 6)

# ╔═╡ 9cadf5aa-bb63-4179-8a99-04714eb74617
substitute(cm(Val(:BiNormal), 6), μ => mean_binormal) |> simplify |> terms

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
Distributions = "~0.25.123"
Latexify = "~0.16.10"
PlutoUI = "~0.7.79"
Symbolics = "~7.14.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.5"
manifest_format = "2.0"
project_hash = "a71ab511810f3da7c2d8f768e19480635cc0b2a1"

[[deps.ADTypes]]
git-tree-sha1 = "f7304359109c768cf32dc5fa2d371565bb63b68a"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.21.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "856ecd7cebb68e5fc87abecd2326ad59f0f911f3"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.43"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "7e35fca2bdfba44d797c53dfe63a51fabf39bfc0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.4.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "d81ae5489e13bc03567d4fbbb06c546a5e53c857"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.22.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = ["CUDSS", "CUDA"]
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bijections]]
git-tree-sha1 = "a2d308fcd4c2fb90e943cf9cd2fbfa9c32b69733"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.2.2"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "fbcc7610f6d8348428f722ecbe0e6cfe22e672c6"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.123"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "c249d86e97a7e8398ce2068dce4c078a1c3464de"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.16"

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"
    DomainSetsRandomExt = "Random"

    [deps.DomainSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "Test"]
git-tree-sha1 = "3f50fa86c968fc1a9e006c07b6bc40ccbb1b704d"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.4"

[[deps.EnumX]]
git-tree-sha1 = "7bebc8aad6ee6217c78c5ddcf7ed289d65d0263e"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.6"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "2f979084d1e13948a3352cf64a25df6bd3b4dca3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.16.0"
weakdeps = ["PDMats", "SparseArrays", "StaticArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IntervalSets]]
git-tree-sha1 = "d966f85b3b7a8e49d034d27a189e9a4874b4391a"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.13"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "53f817d3e84537d84545e0ad749e483412dd6b2a"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "d38b8653b1cdfac5a7da3b819c0a8d6024f9a18c"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.13"

    [deps.MultivariatePolynomials.extensions]
    MultivariatePolynomialsChainRulesCoreExt = "ChainRulesCore"

    [deps.MultivariatePolynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "22df8573f8e7c593ac205455ca088989d0a2c7a0"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.7"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e4cff168707d441cd6bf3ff7e4832bdf34278e4a"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.37"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3ac7038a98ef6977d44adeadc73cc6f596c08109"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.79"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "522f093a29b31a93e34eaea17ba055d850edea28"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.1"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.ReadOnlyArrays]]
git-tree-sha1 = "e6f7ddf48cf141cb312b078ca21cb2d29d0dc11d"
uuid = "988b38a3-91fc-5605-94a2-ee2116b3bd83"
version = "0.2.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "7257165d5477fd1025f7cb656019dcb6b0512c38"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.17"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SciMLPublic]]
git-tree-sha1 = "0ba076dbdce87ba230fff48ca9bca62e1f345c9b"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.0.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5acc6a41b3082920f79ca3c759acbcecf18a8d78"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.7.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eee1b9ad8b29ef0d936e3ec9838c7ec089620308"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.16"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "aceda6f4e598d331548e04cc6b2124a6148138e3"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.10"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91f091a8716a6bb38417a6e6f274602a19aaa685"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.2"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "94c58884e013efff548002e8dc2fdd1cb74dfce5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.46"

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

    [deps.SymbolicIndexingInterface.weakdeps]
    PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils", "TermInterface"]
git-tree-sha1 = "5085671d2cba1eb02136a3d6661c583e801984c1"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "1.1.0"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "EnumX", "ExproniconLite", "LinearAlgebra", "MacroTools", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "PrecompileTools", "ReadOnlyArrays", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "WeakCacheSets"]
git-tree-sha1 = "8a467146b35b676b42b9865567216d5e44f4d504"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "4.18.1"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsChainRulesCoreExt = "ChainRulesCore"
    SymbolicUtilsDistributionsExt = "Distributions"
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "AbstractPlutoDingetjes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "PrecompileTools", "Preferences", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLPublic", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "cf4904e6d73a14d7a662801dd638645aa6efb213"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "7.14.0"

    [deps.Symbolics.extensions]
    SymbolicsD3TreesExt = "D3Trees"
    SymbolicsDistributionsExt = "Distributions"
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsHypergeometricFunctionsExt = "HypergeometricFunctions"
    SymbolicsLatexifyExt = ["Latexify", "LaTeXStrings"]
    SymbolicsLuxExt = "Lux"
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"
    SymbolicsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Symbolics.weakdeps]
    D3Trees = "e3df1716-f71e-5df9-9e2d-98e193103c45"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    HypergeometricFunctions = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    Lux = "b2108857-7c20-44ae-9111-449ecde12c47"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TaskLocalValues]]
git-tree-sha1 = "67e469338d9ce74fc578f7db1736a74d93a49eb8"
uuid = "ed4db957-447d-4319-bfb6-7fa9ae7ecf34"
version = "0.1.3"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.WeakCacheSets]]
git-tree-sha1 = "386050ae4353310d8ff9c228f83b1affca2f7f38"
uuid = "d30d5f5c-d141-4870-aa07-aabb0f5fe7d5"
version = "0.1.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"
"""

# ╔═╡ Cell order:
# ╟─589dae68-0c2d-4c5a-a58d-2d09f0f51c74
# ╟─857d8253-936d-4226-b9af-87553e65e0fb
# ╟─c7a48d7b-4bf7-429e-89de-18bbf5a8a7fd
# ╟─e8d2ef67-e395-4915-9b19-4df3e10327dc
# ╠═a877d102-099f-11f0-37d1-39f0ef1f87e3
# ╠═e5e0e4e2-2df1-4536-9cc5-bdcec6fc13de
# ╠═a357df88-fae3-40ad-b636-a3fbfd5358bb
# ╠═c875992e-a423-4376-9c58-88fc7a7c1a8c
# ╟─b6466bc4-ff63-4d88-9913-c2b1b80ea9eb
# ╠═1861682e-b53c-4acb-88b6-7a20eaa5a43c
# ╟─f1c4bf46-c020-4bae-a805-0e834c90c82d
# ╟─d5fcbe3c-ff12-4394-a225-16d211024693
# ╠═2d52e729-1c62-49a3-8acc-a9e19a3c2797
# ╠═68746676-7121-4ef4-afbf-ded292130da6
# ╟─30de8ac9-3935-4e0a-b287-94884d8d30e8
# ╟─afc8be14-3e99-4eab-bbbd-673ec87e1052
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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
