# Transform Shiller's ie_data.xlsx into the annual ISSM building blocks (`load_shiller_annual`)
# and save them as CSV.
#
# Usage (from the repo root):
#   julia --project analisys/Ibbotson-Sinquefield-SimulationModel/export_shiller_annual.jl [input.xlsx] [output.csv]

using CSV
using TimeSeries: colnames
using MacroFinanceScenarios.IbbotsonSinquefield

const ROOT = joinpath(@__DIR__, "..", "..")

input  = length(ARGS) >= 1 ? ARGS[1] : joinpath(ROOT, "data", "ie_data.xlsx")
output = length(ARGS) >= 2 ? ARGS[2] : joinpath(ROOT, "data", "shiller_annual.csv")

ta = load_shiller_annual(input)
CSV.write(output, ta)
println("Wrote $(length(ta)) rows × $(length(colnames(ta))) columns to $output")
