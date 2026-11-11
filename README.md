# TwoScalesVolatility.jl

Two-Scales Realized Volatility.

The two-scales realized volatility of Zhang, Mykland and Aït-Sahalia (2005): the integrated variance of a day from prices observed every second with market microstructure noise, the average of the realized variances of K sub-grids less a multiple of the realized variance of every tick, with the small-sample adjustment; the noise variance; the optimal sparse step, sub-grid size and K; and the Monte Carlo study of their Table 2 on the Heston model.

Article: Zhang, L., Mykland, P. A. and Aït-Sahalia, Y. (2005). A tale of two time scales: determining integrated volatility with noisy high-frequency data. Journal of the American Statistical Association 100(472), 1394–1411. [doi:10.1198/016214505000000169](https://doi.org/10.1198/016214505000000169)

## Usage

```julia
using TwoScalesVolatility, Random

# A day of 23,400 one-second log prices with a daily volatility of 1% and noise of 5 basis points; with tick data,
# the log prices of a Tables.jl table instead: y = log.(ticks.price)
y = simulate_constant_day(Xoshiro(1); volatility = 0.01, noise = 0.0005)
n = length(y) - 1

sqrt(realized_variance(y))      # every tick, the fifth-best estimator: about 0.108, ten times the true 0.01
sqrt(noise_variance(y))         # the noise, equation (22): about 0.0005

# K = c* n^(2/3), equation (10), from the noise variance and a first guess at the integrated variance from five-minute
# returns (a day is T = 1, and a constant volatility has a quarticity of its integrated variance squared)
first_guess = sparse_variance(y, 300)
K = round(Int, optimal_c(1.0, noise_variance(y), first_guess^2) * n^(2 / 3))
sqrt(two_scales_variance(y, K)) # the two-scales estimator, the first-best: about 0.0098
sqrt(noise_variance(y, K))      # the noise, equation (65): about 0.0005

# Table 2 of the article, a 9×5 matrix, on 1,000 simulated days of the Heston model of Section 7.1 (the article's 25,000
# is the default)
table2(; paths = 1_000)
```

## Vignette

[Two-Scales Realized Volatility](vignette/TwoScalesVolatility.tex), its LaTeX source: Pollis compiles it in Help > Documentation > Package Vignettes.

## License

MIT (see [LICENSE](LICENSE)).

## Tests

The tests compute the published tables below from their data and check every printed value, rounded to the digits printed:

- Table 2, page 1405 (typed in by a coding agent) (`test/data/targets/table2.csv`)

They also run the tests of `test/tests/tests.jl`.

Run them with `using Pkg; Pkg.test("TwoScalesVolatility")`.

Made with [Pollis](https://pollis.org).
