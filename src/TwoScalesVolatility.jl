module TwoScalesVolatility

using Distributions, Random, Statistics

export HestonDesign, simulate_day, realized_variance, sparse_variance, average_variance, two_scales_variance,
	noise_variance, optimal_sparse, optimal_average, optimal_c, table2, monte_carlo_table2, simulate_constant_day,
	recovery

"""
    HestonDesign(; mu = 0.05, kappa = 5.0, alpha = 0.04, gamma = 0.5, rho = -0.5, noise = 0.0005, T = 1 / 252, n = 23_400)

The design of the Monte Carlo study of Section 7.1 of Zhang, Mykland and Aït-Sahalia (2005): the Heston model
dX = (mu - v/2) dt + sqrt(v) dB, dv = kappa (alpha - v) dt + gamma sqrt(v) dW, corr(B, W) = rho, annualized; a day
of `T` years observed `n` times at equal steps (one second over 6.5 hours), each observation of the log price
carrying independent Gaussian noise of standard deviation `noise`.
"""
Base.@kwdef struct HestonDesign
	mu::Float64 = 0.05
	kappa::Float64 = 5.0
	alpha::Float64 = 0.04
	gamma::Float64 = 0.5
	rho::Float64 = -0.5
	noise::Float64 = 0.0005
	T::Float64 = 1 / 252
	n::Int = 23_400
end

"""
    simulate_day(design, rng; x = Vector{Float64}(undef, design.n + 1), y = similar(x))

One day of the design by the Euler scheme at its step (variance truncated at zero), the variance starting from its
stationary Gamma distribution (shape 2 kappa alpha / gamma^2, scale gamma^2 / (2 kappa)): the article does not say how
it starts, and this start is the one its asymptotic rows fit. Fills `x`, the efficient log prices at the n + 1 times,
and `y`, the observed ones, `x` plus the noise, and returns `(; iv, iq)`, the integrated variance and quarticity of
the day, the integrals of v and v^2, as Riemann sums on the grid.
"""
function simulate_day(d::HestonDesign, rng::AbstractRNG; x = Vector{Float64}(undef, d.n + 1), y = similar(x))
	dt = d.T / d.n
	v = rand(rng, Gamma(2 * d.kappa * d.alpha / d.gamma^2, d.gamma^2 / (2 * d.kappa)))
	orthogonal = sqrt(1 - d.rho^2)
	x[1] = 0.0
	iv = 0.0
	iq = 0.0
	for i in 1:d.n
		positive = max(v, 0.0)
		iv += positive * dt
		iq += positive^2 * dt
		z1 = randn(rng)
		z2 = randn(rng)
		root = sqrt(positive * dt)
		x[i + 1] = x[i] + (d.mu - positive / 2) * dt + root * z1
		v += d.kappa * (d.alpha - positive) * dt + d.gamma * root * (d.rho * z1 + orthogonal * z2)
	end
	for i in eachindex(x)
		y[i] = x[i] + d.noise * randn(rng)
	end
	return (; iv, iq)
end

"""
    realized_variance(y)

The realized variance of the log prices `y`, the sum of their squared returns: [Y, Y]^(all) of the article, its
fifth-best estimator.
"""
realized_variance(y::AbstractVector) = sum(abs2(y[i + 1] - y[i]) for i in 1:length(y) - 1)

"""
    sparse_variance(y, step; whole_day = false)

The realized variance of every `step`-th log price of `y`, from the first: [Y, Y]^(sparse) of the article, its
fourth-best estimator (third-best with the optimal step), with `length(y) - 1` divided by `step` returns. When `step`
does not divide the n returns, the last ticks of the day are left out; with `whole_day`, the sum is multiplied by n
over the ticks its returns span, so that it covers the whole day, as the small-sample biases of Table 2 show the
article's estimator does (the article does not say).
"""
function sparse_variance(y::AbstractVector, step::Integer; whole_day::Bool = false)
	n = length(y) - 1
	returns = n ÷ step
	total = sum(abs2(y[i + step] - y[i]) for i in 1:step:step * (returns - 1) + 1)
	return whole_day ? total * n / (returns * step) : total
end

"""
    average_variance(y, K; whole_day = false)

The average over the `K` regular sub-grids of `y` (the k-th holds the log prices k, k + K, k + 2K, ...) of their
realized variances: [Y, Y]^(avg) of the article, equation (32), its second-best estimator. The returns of all the
sub-grids are the K-step returns of `y`, so the average is their sum of squares divided by K. Each sub-grid leaves out
some ticks at the start or the end of the day, K - 1 on average; with `whole_day`, the average is multiplied by
n / (n - K + 1), so that it covers the whole day, as the small-sample biases of Table 2 show the article's estimator
does (the article does not say).
"""
function average_variance(y::AbstractVector, K::Integer; whole_day::Bool = false)
	n = length(y) - 1
	total = sum(abs2(y[i + K] - y[i]) for i in 1:n + 1 - K) / K
	return whole_day ? total * n / (n - K + 1) : total
end

"""
    two_scales_variance(y, K; adjust = true)

The two-scales estimator of the integrated variance from the log prices `y` with `K` sub-grids, the first-best of
the article: [Y, Y]^(avg) - (n̄ / n) [Y, Y]^(all), equation (55), with n̄ = (n - K + 1) / K, equation (34); with
`adjust`, divided by 1 - n̄ / n, the small-sample adjustment of equation (64).
"""
function two_scales_variance(y::AbstractVector, K::Integer; adjust::Bool = true)
	n = length(y) - 1
	nbar = (n - K + 1) / K
	estimate = average_variance(y, K) - nbar / n * realized_variance(y)
	return adjust ? estimate / (1 - nbar / n) : estimate
end

"""
    noise_variance(y)

The estimate of the variance of the noise from all the log prices `y`: [Y, Y]^(all) / (2n), equation (22).
"""
noise_variance(y::AbstractVector) = realized_variance(y) / (2 * (length(y) - 1))

"""
    noise_variance(y, K)

The adjusted estimate of the variance of the noise from the log prices `y` and `K` sub-grids,
([Y, Y]^(all) - [Y, Y]^(avg)) / (2 (n - n̄)), equation (65), unbiased to a higher order than [Y, Y]^(all) / (2n).
"""
function noise_variance(y::AbstractVector, K::Integer)
	n = length(y) - 1
	nbar = (n - K + 1) / K
	return (realized_variance(y) - average_variance(y, K)) / (2 * (n - nbar))
end

"""
    optimal_sparse(T, noise_variance, quarticity)

The optimal number of sparse returns of the third-best estimator for equidistant sampling,
n*_sparse = (T / (4 (E eps^2)^2) * quarticity)^(1/3) (Section 1.2.3), `quarticity` the integral of sigma^4.
"""
optimal_sparse(T, noise_variance, quarticity) = cbrt(T / (4 * noise_variance^2) * quarticity)

"""
    optimal_average(T, noise_variance, quarticity)

The optimal average size of a sub-grid of the second-best estimator, n̄* = (T / (6 (E eps^2)^2) * quarticity)^(1/3),
equation (8); the number of sub-grids is then K* = n / n̄*.
"""
optimal_average(T, noise_variance, quarticity) = cbrt(T / (6 * noise_variance^2) * quarticity)

"""
    optimal_c(T, noise_variance, quarticity)

The optimal constant of the first-best estimator, c* = (T / (12 (E eps^2)^2) * quarticity)^(-1/3), equation (10);
the number of sub-grids is then K = c* n^(2/3).
"""
optimal_c(T, noise_variance, quarticity) = cbrt(T / (12 * noise_variance^2) * quarticity)^-1

"""
The variance of the errors within bins of the paths by integrated variance, then quarticity, averaged with the bins'
sizes, as for (73), and its standard error, from the fourth central moment within each bin.
"""
function binned_variance(errors, iv, iq, bins)
	first, second = bins
	order = sortperm(iv)
	total = 0.0
	spread = 0.0
	count = 0
	for group in Iterators.partition(order, cld(length(order), first))
		inner = group[sortperm(iq[group])]
		for bin in Iterators.partition(inner, cld(length(inner), second))
			length(bin) < 2 && continue
			values = errors[bin]
			s2 = var(values)
			m4 = mean(abs2.(abs2.(values .- mean(values))))
			total += s2 * length(bin)
			spread += length(bin) * max(m4 - s2^2, 0.0)
			count += length(bin)
		end
	end
	return total / count, sqrt(spread) / count
end

"The mean of `values` and its standard error."
mean_se(values) = (mean(values), std(values) / sqrt(length(values)))

"The square root of `bias`^2 + `variance` and its standard error, by the delta method, the two taken as independent."
function rmse_se((bias, bias_se), (variance, variance_se))
	value = sqrt(bias^2 + variance)
	return value, sqrt((bias * bias_se)^2 + (variance_se / 2)^2) / value
end

"""
    simulate_constant_day(rng; volatility = 0.01, noise = 0.0, n = 23_400, x = Vector{Float64}(undef, n + 1), y = similar(x))

A day of log prices with a known constant volatility: `n` returns of a Brownian motion without drift whose
integrated variance over the day is `volatility`^2, observed with independent Gaussian noise of standard deviation
`noise`. Fills `x`, the efficient log prices at the n + 1 times, and `y`, the observed ones, and returns `y`.
"""
function simulate_constant_day(rng::AbstractRNG; volatility::Real = 0.01, noise::Real = 0.0, n::Integer = 23_400,
		x = Vector{Float64}(undef, n + 1), y = similar(x))
	step = volatility / sqrt(n)
	x[1] = 0.0
	for i in 1:n
		x[i + 1] = x[i] + step * randn(rng)
	end
	for i in eachindex(x)
		y[i] = x[i] + noise * randn(rng)
	end
	return y
end

"""
    recovery(; volatility = 0.01, noises = [0.0, 0.0005, 0.001, 0.002], days = 400, n = 23_400, seed = 7)

How well the two-scales estimator recovers a known integrated variance, and equation (65) a known noise variance:
for each noise standard deviation of `noises` (0 to 20 basis points), `days` days of `simulate_constant_day` with
daily `volatility`, the first-best estimator with K = c* n^(2/3) from the true quarticity and noise variance, at
least 2. Returns a vector of named tuples, one for each noise: `noise`; `K`; `variance`, the true integrated variance
volatility^2; `expected`, what the estimator gives on average, variance × (1 - K / (n + 1)), since its sub-grids leave
out ticks at the ends of the day; `estimate` and `estimate_se`, the mean of the estimates over the days and its
standard error; `noise_variance`, the true noise variance; and `noise_estimate` and `noise_se`, the mean of the
estimates of equation (65) and its standard error.
"""
function recovery(; volatility::Real = 0.01, noises = [0.0, 0.0005, 0.001, 0.002], days::Integer = 400, n::Integer = 23_400, seed::Integer = 7)
	variance = volatility^2
	buffers = [(x = Vector{Float64}(undef, n + 1), y = Vector{Float64}(undef, n + 1)) for _ in 1:Threads.maxthreadid()]
	return map(enumerate(noises)) do (level, noise)
		# With T = 1, the quarticity of a constant volatility is its integrated variance squared
		K = noise == 0 ? 2 : clamp(round(Int, optimal_c(1.0, noise^2, variance^2) * n^(2 / 3)), 2, n ÷ 2)
		estimates = zeros(days)
		noise_estimates = zeros(days)
		Threads.@threads :static for day in 1:days
			x, y = buffers[Threads.threadid()]
			simulate_constant_day(Xoshiro(UInt64[seed, level, day]); volatility, noise, n, x, y)
			estimates[day] = two_scales_variance(y, K)
			noise_estimates[day] = noise_variance(y, K)
		end
		estimate, estimate_se = mean_se(estimates)
		noise_estimate, noise_se = mean_se(noise_estimates)
		(; noise, K, variance, expected = variance * (1 - K / (n + 1)), estimate, estimate_se,
			noise_variance = noise^2, noise_estimate, noise_se)
	end
end

"""
    table2(design = HestonDesign(); paths = 25_000, seed = 2005, sparse_step = 300, bins = (50, 5))

Table 2 of Zhang, Mykland and Aït-Sahalia (2005): the Monte Carlo study of the five estimators of the integrated
variance, simulated with `simulate_day` on `paths` days. Returns a 9×5 matrix: the rows small-sample bias,
asymptotic bias, small-sample variance, asymptotic variance, small-sample RMSE, asymptotic RMSE, small-sample
relative bias, relative variance and relative RMSE; the columns the fifth-best (every tick), fourth-best (every
`sparse_step` ticks, five minutes), third-best (the optimal sparse step), second-best (sub-grids of the optimal
average size) and first-best (two scales, adjusted, K = c* n^(2/3)) estimators. The optimal choices of each path use
its integrated quarticity and the noise variance of the design, a sparse step and a number of sub-grids kept
within the day.

The small-sample bias is the mean error over the paths, its variance the variance of the errors within `bins` of the
paths by integrated variance, then quarticity, averaged across the bins, as the article says without giving the bins;
the RMSE is the square root of bias squared plus variance; the relative rows are the same for the errors divided by
the integrated variance. The asymptotic rows average over the paths the bias and the variance the theory predicts for
each, its noise part with the second-order terms of equations (20), (37) and (60), its RMSE the square root of their
mean bias squared plus mean variance. The sparse and averaged estimators cover the whole day (`whole_day`), with
n / step sparse returns and sub-grids of n / K returns in the theory's bias and variance. Each path has its own
random numbers, from `seed` and its number, so the result does not depend on the number of threads.
"""
table2(d::HestonDesign = HestonDesign(); options...) = monte_carlo_table2(d; options...).values

"""
    monte_carlo_table2(design = HestonDesign(); paths = 25_000, seed = 2005, sparse_step = 300, bins = (50, 5))

Table 2 as `table2` computes it, with the Monte Carlo standard error of each value: returns `(; values,
standard_errors)`, two 9×5 matrices. A mean over the paths has the standard deviation over the paths divided by
sqrt(paths); a binned variance the standard error of the variance within each bin, from its fourth central moment,
combined across the bins; an RMSE the standard error from those of its bias and variance by the delta method. Two runs
of the same size, the article's and this one, differ by about sqrt(2) standard errors.
"""
function monte_carlo_table2(d::HestonDesign = HestonDesign(); paths::Integer = 25_000, seed::Integer = 2005, sparse_step::Integer = 300, bins = (50, 5))
	n = d.n
	noise2 = d.noise^2
	noise4 = 3 * noise2^2
	noisevar2 = 2 * noise2^2
	estimators = 5
	errors = zeros(paths, estimators)
	bias = zeros(paths, estimators)
	variance = zeros(paths, estimators)
	iv = zeros(paths)
	iq = zeros(paths)
	buffers = [(x = Vector{Float64}(undef, n + 1), y = Vector{Float64}(undef, n + 1)) for _ in 1:Threads.maxthreadid()]
	Threads.@threads :static for m in 1:paths
		x, y = buffers[Threads.threadid()]
		day = simulate_day(d, Xoshiro(UInt64[seed, m]); x, y)
		iv[m] = day.iv
		iq[m] = day.iq
		discretization = d.T * day.iq
		second(xx) = 8 * xx * noise2 - 2 * noisevar2
		# Fifth-best: every tick
		errors[m, 1] = realized_variance(y) - day.iv
		bias[m, 1] = 2 * n * noise2
		variance[m, 1] = 4 * n * noise4 + second(realized_variance(x)) + 2 * discretization / n
		# Fourth- and third-best: sparse, every sparse_step ticks or at the optimal step
		optimal = clamp(round(Int, n / optimal_sparse(d.T, noise2, day.iq)), 1, n)
		for (column, step) in ((2, sparse_step), (3, optimal))
			nsparse = n / step
			errors[m, column] = sparse_variance(y, step; whole_day = true) - day.iv
			bias[m, column] = 2 * nsparse * noise2
			variance[m, column] = 4 * nsparse * noise4 + second(sparse_variance(x, step; whole_day = true)) + 2 * discretization / nsparse
		end
		# Second-best: the average over K sub-grids of the optimal average size
		K = clamp(round(Int, n / optimal_average(d.T, noise2, day.iq)), 1, n)
		nbar = n / K
		errors[m, 4] = average_variance(y, K; whole_day = true) - day.iv
		bias[m, 4] = 2 * nbar * noise2
		variance[m, 4] = nbar / K * 4 * noise4 + second(average_variance(x, K; whole_day = true)) / K + 4 * discretization / (3 * nbar)
		# First-best: two scales, adjusted, K = c* n^(2/3)
		K = clamp(round(Int, optimal_c(d.T, noise2, day.iq) * n^(2 / 3)), 2, n ÷ 2)
		c = K / n^(2 / 3)
		errors[m, 5] = two_scales_variance(y, K) - day.iv
		bias[m, 5] = 0.0
		variance[m, 5] = n^(-1 / 3) * (8 * noise2^2 / c^2 + c * 4 * discretization / 3) + n^(-2 / 3) / c * second(average_variance(x, K))
	end
	relative = errors ./ iv
	values = zeros(9, estimators)
	standard_errors = zeros(9, estimators)
	for j in 1:estimators
		small = mean_se(view(errors, :, j))
		spread = binned_variance(view(errors, :, j), iv, iq, bins)
		asymptotic = mean_se(view(bias, :, j))
		predicted = mean_se(view(variance, :, j))
		relative_bias = mean_se(view(relative, :, j))
		relative_spread = binned_variance(view(relative, :, j), iv, iq, bins)
		column = [small, asymptotic, spread, predicted, rmse_se(small, spread), rmse_se(asymptotic, predicted),
			relative_bias, relative_spread, rmse_se(relative_bias, relative_spread)]
		values[:, j] = first.(column)
		standard_errors[:, j] = last.(column)
	end
	return (; values, standard_errors)
end

end
