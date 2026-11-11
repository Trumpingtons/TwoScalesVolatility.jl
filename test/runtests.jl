# Made by Pollis (Package and Publish, on the page of an article method) from method.toml: change the method, not this file.
using Test, CSV, Downloads, Tables, TwoScalesVolatility

const DATA = joinpath(@__DIR__, "data")

"An input table: a file of the tests, or downloaded from its link into memory."
function input(source; delim = ',', columns = nothing, header = true, skipto = nothing, limit = nothing)
	bytes = startswith(source, r"https?://") ? take!(Downloads.download(source, IOBuffer())) : read(joinpath(DATA, source))
	named = columns === nothing && header
	header = columns !== nothing ? Symbol.(columns) : named ? 1 : false
	rows = (; skipto = something(skipto, named ? 2 : 1), limit)
	file = delim == "whitespace" ? CSV.File(bytes; header, rows..., delim = ' ', ignorerepeated = true, stripwhitespace = true) : CSV.File(bytes; header, rows..., delim = delim == "tab" ? '\t' : only(delim))
	return Tables.columntable(file)
end

"The published table: its row labels and its values as printed, `nothing` where a cell is empty."
function target(file)
	table = CSV.File(joinpath(DATA, "targets", file * ".csv"); types = String, delim = ',')
	columns = Tables.columnnames(table)
	rows = [string(coalesce(row[columns[1]], "")) for row in table]
	values = [let cell = row[column]; cell === missing ? nothing : strip(cell) end for row in table, column in columns[2:end]]
	return rows, string.(columns[2:end]), values
end

"Tests each printed value against the result rounded to the digits printed, but the untested cells, each (row label, column name)."
function check(result, file, tolerance, untested = Tuple{String, String}[])
	rows, columns, values = target(file)
	result = result isa AbstractVecOrMat ? result : Tables.matrix(result)
	result = reshape(result, size(result, 1), :)
	@test size(result) == size(values)
	@test isempty(setdiff(untested, [(row, column) for row in rows for column in columns]))
	for i in axes(values, 1), j in axes(values, 2)
		printed = values[i, j]
		(printed === nothing || isempty(printed) || (rows[i], columns[j]) in untested) && continue
		mantissa, exponent = occursin(r"[eE]", printed) ? split(printed, r"[eE]") : (printed, "0")
		digits = (occursin('.', mantissa) ? length(split(mantissa, '.')[2]) : 0) - parse(Int, exponent)
		published = parse(Float64, printed)
		computed = round(result[i, j]; digits)
		@testset "$(rows[i]), $(columns[j])" begin
			@test isapprox(computed, published; atol = (tolerance + 1e-6) * 10.0^-digits)
		end
	end
end

@testset "TwoScalesVolatility" begin
	include("tests/tests.jl")
	@testset "Table 2, page 1405 (typed in by a coding agent)" begin
		check(table2(), "table2", 1, [("Small-sample bias", "Fifth-best"), ("Small-sample bias", "Fourth-best"), ("Small-sample bias", "Third-best"), ("Small-sample bias", "Second-best"), ("Small-sample bias", "First-best"), ("Asymptotic bias", "Third-best"), ("Asymptotic bias", "Second-best"), ("Small-sample variance", "Fifth-best"), ("Small-sample variance", "Fourth-best"), ("Small-sample variance", "Third-best"), ("Small-sample variance", "Second-best"), ("Small-sample variance", "First-best"), ("Asymptotic variance", "Fifth-best"), ("Asymptotic variance", "Fourth-best"), ("Asymptotic variance", "Third-best"), ("Asymptotic variance", "Second-best"), ("Asymptotic variance", "First-best"), ("Small-sample RMSE", "Fifth-best"), ("Small-sample RMSE", "Fourth-best"), ("Small-sample RMSE", "Third-best"), ("Small-sample RMSE", "Second-best"), ("Small-sample RMSE", "First-best"), ("Asymptotic RMSE", "Fourth-best"), ("Asymptotic RMSE", "Third-best"), ("Asymptotic RMSE", "Second-best"), ("Asymptotic RMSE", "First-best"), ("Small-sample relative bias", "Fifth-best"), ("Small-sample relative bias", "Fourth-best"), ("Small-sample relative bias", "Third-best"), ("Small-sample relative bias", "Second-best"), ("Small-sample relative bias", "First-best"), ("Small-sample relative variance", "Fifth-best"), ("Small-sample relative variance", "Fourth-best"), ("Small-sample relative variance", "Third-best"), ("Small-sample relative variance", "Second-best"), ("Small-sample relative variance", "First-best"), ("Small-sample relative RMSE", "Fifth-best"), ("Small-sample relative RMSE", "Fourth-best"), ("Small-sample relative RMSE", "Third-best"), ("Small-sample relative RMSE", "Second-best"), ("Small-sample relative RMSE", "First-best")])
	end
end
