# The two-scales estimator against rTSCov of the R package highfrequency, an independent implementation of equation
# (64) with J = 1: on days with a known volatility and noise of 5 to 20 basis points, the two agree within 1e-5 of
# their value. They differ by a few millionths because highfrequency counts the n + 1 prices where the article
# counts the n returns, in n̄ = (n - K + 1) / K and in the n of n̄ / n. Optional: it runs only when named,
# Pkg.test("TwoScalesVolatility"; test_args = ["rtscov"]), and is skipped when R, highfrequency or xts is missing.

using Random

"Whether Rscript runs and has highfrequency and xts."
function has_highfrequency()
	Sys.which("Rscript") === nothing && return false
	return success(pipeline(`Rscript -e "quit(status = as.integer(!all(sapply(c('highfrequency', 'xts'), requireNamespace, quietly = TRUE))))"`; stdout = devnull, stderr = devnull))
end

"""The rTSCov test: on four days, the two-scales estimator against rTSCov, within 1e-5 of their value."""
function test_rtscov()
	@testset "rTSCov of R highfrequency" begin
		if !has_highfrequency()
			@info "Not tested: rTSCov of R highfrequency, since Rscript, highfrequency or xts is missing."
			return
		end
		mktempdir() do folder
			for (case, (noise, K)) in enumerate([(0.0005, 35), (0.001, 87), (0.002, 219), (0.0005, 300)])
				y = simulate_constant_day(Xoshiro(case); volatility = 0.01, noise)
				file = joinpath(folder, "day$case.txt")
				write(file, join(string.(exp.(log(100) .+ y)), "\n"))
				script = "suppressMessages({library(highfrequency); library(xts)}); p <- scan('$(escape_string(file))', quiet = TRUE); " *
					"x <- xts(p, as.POSIXct('2005-01-03 09:30:00', tz = 'UTC') + seq_along(p) - 1); cat(sprintf('%.17g', rTSCov(x, K = $K, J = 1)))"
				r = parse(Float64, readchomp(pipeline(`Rscript -e $script`; stderr = devnull)))
				ours = two_scales_variance(y, K)
				@testset "Noise $(round(noise * 1e4; digits = 1)) bp, K = $K: $(round(ours; sigdigits = 8)) against rTSCov $(round(r; sigdigits = 8))" begin
					@test abs(ours / r - 1) < 1e-5
				end
			end
		end
	end
end

if "rtscov" in ARGS
	test_rtscov()
else
	@info "Not tested: rTSCov of R highfrequency, which needs R. To test it: Pkg.test(\"TwoScalesVolatility\"; test_args = [\"rtscov\"])"
end
