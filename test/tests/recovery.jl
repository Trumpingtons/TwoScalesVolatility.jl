# Recovery: days with a known constant volatility of 1% a day and noise of 0 to 20 basis points. The mean two-scales
# estimate must be within 3 standard errors of what the estimator gives on average, the integrated variance times
# 1 - K / (n + 1) (its sub-grids leave out ticks at the ends of the day), and within 2% of the integrated variance
# itself; the mean estimate of the noise variance, equation (65), within 3 standard errors of the true one.

@testset "Recovery of a known integrated variance and noise variance" begin
	for r in recovery()
		@testset "Noise $(round(r.noise * 1e4; digits = 1)) bp, K = $(r.K): integrated variance $(round(r.estimate; sigdigits = 5)) against $(r.variance), noise variance $(round(r.noise_estimate; sigdigits = 4)) against $(r.noise_variance)" begin
			@test abs(r.estimate - r.expected) < 3 * r.estimate_se
			@test abs(r.estimate / r.variance - 1) < 0.02
			@test abs(r.noise_estimate - r.noise_variance) < 3 * r.noise_se
		end
	end
end
