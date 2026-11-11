# Table 2 of Zhang, Mykland and Aït-Sahalia (2005) within the noise of the simulation. The article's values come from
# one run of 25,000 simulated days, these from another, so neither gives the other's digits. Each value is compared as
# Pollis.compare compares random values: the gap beyond the rounding of the printed value, half a unit of its last
# printed digit, in standard errors of the difference of two runs of the same size, sqrt(2) times this run's standard
# error, must be under 3: Pollis.compare calls a value different from 2, which about 2 of 40 correct values pass by
# chance, and 3 about 1 in 370. The values the theory gives without simulation are tested to the digits printed in the
# target table, the fifth-best asymptotic RMSE among them; the fifth-best relative variance and relative RMSE are not tested: E[1 / v^2] is infinite for the
# Gamma distribution of v (shape 1.6), so their sample values never settle.

const TABLE2_EXACT = [("Asymptotic bias", "Fifth-best"), ("Asymptotic bias", "Fourth-best"), ("Asymptotic bias", "First-best"), ("Asymptotic RMSE", "Fifth-best")]
const TABLE2_UNSTABLE = [("Small-sample relative variance", "Fifth-best"), ("Small-sample relative RMSE", "Fifth-best")]

@testset "Table 2, page 1405, within the noise of the simulation" begin
	(; values, standard_errors) = monte_carlo_table2()
	rows, columns, printed = target("table2")
	for i in axes(printed, 1), j in axes(printed, 2)
		cell = (rows[i], columns[j])
		(cell in TABLE2_EXACT || cell in TABLE2_UNSTABLE) && continue
		mantissa, exponent = occursin(r"[eE]", printed[i, j]) ? split(printed[i, j], r"[eE]") : (printed[i, j], "0")
		decimals = (occursin('.', mantissa) ? length(split(mantissa, '.')[2]) : 0) - parse(Int, exponent)
		gap = values[i, j] - parse(Float64, printed[i, j])
		z = max(abs(gap) - 0.5 * 10.0^-decimals, 0.0) / (sqrt(2) * standard_errors[i, j])
		@testset "$(rows[i]), $(columns[j]): $(round(values[i, j]; sigdigits = 5)) against $(printed[i, j]), z = $(round(z; digits = 2))" begin
			@test z < 3
		end
	end
end
