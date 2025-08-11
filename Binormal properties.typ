#import "@preview/physica:0.9.4": *

= Bi-normal distribution

== Probability distribution function
$
    f(x; λ, μ_1, σ_1, μ_2, σ_2)
    = λ N(x; μ_1, σ_1) + (1 - λ) N(x; μ_2, σ_2)
$
where $display(N(x; μ, σ) = 1/(sqrt(2pi)σ) exp(-(x-μ)^2/(2 σ^2)))$

== Cumulative distribution function
$
    F(x; λ, μ_1, σ_1, μ_2, σ_2)
    = λ Phi (x; μ_1, σ_1) + (1 - λ) Phi (x; μ_2, σ_2)
$
where
$
    Phi (x; μ, σ) = 1/2 [1 + erf((x - μ)/(sqrt(2)σ))]
$
is the cdf of the normal distribution and $erf(x) = 2/sqrt(pi) integral_0^x e^(-t^2) d t$.

== Median
The median $m$ satisfies any of the following equivalent equations:
$
  F(m) = 1/2 \
  λ Phi (m; μ_1, σ_1) + (1 - λ) Phi (m; μ_2, σ_2) = 1/2 \
  λ erf((m - μ_1)/σ_1) + (1 - λ) erf((m - μ_2) / σ_2) = 0 \
  erf((m - μ_1)\/σ_1)/erf((m - μ_2)\/σ_2) = 1 - 1/λ
$
