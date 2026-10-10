# Estimated specification on disk: the model parameters as YAML (`save_issm_params`) and
# residual matrices as CSV (`save_residuals`), everything in decimal log units. The EDA
# notebook writes them, the report reads them back without refitting.

const ISSM_NAMES = [:π, :r, :tp, :erp]
const AR1_FIELDS = (:n, :c, :se_c, :phi, :se_phi, :sigma, :r2)

yaml_value(x::Real) = repr(x isa Integer ? x : Float64(x))   # shortest round-trip repr
yaml_list(xs) = "[" * join(yaml_value.(xs), ", ") * "]"

"""
    save_issm_params(path; ar_π, ar_r, ar_tp, erp_target, boot_period, y0, last_year)

Write the Ibbotson–Sinquefield specification to the YAML file `path`, in decimals: the
estimation windows, `erp_target`, the starting point `y0` (state `[π, r, tp, erp]` in
`last_year`), the VAR intercepts `c` and lag diagonal `phi` from `is_coefficients`, and the
AR(1) fits (`n, c, se_c, phi, se_phi, sigma, r2`) for presentation. The AR window is read
from the fits, which must share it. `load_issm_params` reads the file back.
"""
function save_issm_params(path::AbstractString; ar_π::AR1, ar_r::AR1, ar_tp::AR1,
                          erp_target::Real, boot_period::Tuple{Integer,Integer},
                          y0::AbstractVector{<:Real}, last_year::Integer)
    fits = (π = ar_π, r = ar_r, tp = ar_tp)
    allequal(f.period for f in fits) ||
        throw(ArgumentError("the AR(1) fits have different periods: $([f.period for f in fits])"))
    length(y0) == length(ISSM_NAMES) ||
        throw(DimensionMismatch("y0 has length $(length(y0)), expected $(length(ISSM_NAMES))"))
    c, A, _ = is_coefficients(ar_π, ar_r, ar_tp; erp_target)

    open(path, "w") do io
        println(io, "# Ibbotson–Sinquefield specification, decimal log units (save_issm_params)")
        println(io, "ar_period: ", yaml_list(ar_π.period))
        println(io, "boot_period: ", yaml_list(boot_period))
        println(io, "erp_target: ", yaml_value(erp_target))
        println(io, "last_year: ", yaml_value(last_year), "  # simulation starts in last_year + 1")
        println(io, "variables: [", join(ISSM_NAMES, ", "), "]")
        println(io, "y0: ", yaml_list(y0), "  # state in last_year")
        println(io, "c: ", yaml_list(c), "  # VAR intercepts")
        println(io, "phi: ", yaml_list(diag(A)), "  # diagonal of A")
        println(io, "ar1:  # AR(1) fits, for presentation")
        for (name, f) in pairs(fits)
            vals = (length(f.resid), f.c, f.se[1], f.φ, f.se[2], f.σ, f.r2)
            println(io, "  ", name, ": {",
                    join(("$k: $(yaml_value(v))" for (k, v) in zip(AR1_FIELDS, vals)), ", "), "}")
        end
    end
    return path
end

"""
    load_issm_params(path) -> NamedTuple

Read a specification written by `save_issm_params`. Fields: `names` (state variables),
`c`, `A` (diagonal lag matrix), `y0`, `start_year = last_year + 1`, `last_year`,
`ar_period`, `boot_period`, `erp_target` and `ar1`, a vector of NamedTuples
`(variable, n, c, se_c, phi, se_phi, sigma, r2)`, one per AR(1) (a Tables.jl table).
"""
function load_issm_params(path::AbstractString)
    d = YAML.load_file(path)
    names = Symbol.(d["variables"])
    c = Float64.(d["c"])
    A = Matrix{Float64}(Diagonal(Float64.(d["phi"])))
    check_var_dims(c, A, names)
    period(x) = (Int(x[1]), Int(x[2]))
    ar1 = [(; variable = String(v), (k => d["ar1"][v][String(k)] for k in AR1_FIELDS)...)
           for v in ("π", "r", "tp")]
    return (; names, c, A, y0 = Float64.(d["y0"]), start_year = Int(d["last_year"]) + 1,
            last_year = Int(d["last_year"]), ar_period = period(d["ar_period"]),
            boot_period = period(d["boot_period"]), erp_target = Float64(d["erp_target"]), ar1)
end

"""
    save_residuals(path, years, M, names)

Write the residual matrix `M` (one row per year in `years`, columns `names`) to the CSV file
`path`, with a leading `year` column. `load_residuals` reads it back.
"""
function save_residuals(path::AbstractString, years::AbstractVector{<:Integer}, M::AbstractMatrix{<:Real},
                        names::AbstractVector{Symbol})
    size(M) == (length(years), length(names)) ||
        throw(DimensionMismatch("M is $(size(M)), expected $((length(years), length(names)))"))
    cols = (; :year => collect(years), (n => M[:, j] for (j, n) in enumerate(names))...)
    CSV.write(path, cols)
    return path
end

"""
    load_residuals(path) -> (years, M, names)

Read a residual matrix written by `save_residuals`: the years, the `n_years × n` matrix and
the column names.
"""
function load_residuals(path::AbstractString)
    f = CSV.File(path)
    names = filter(!=(:year), propertynames(f))
    M = stack(Vector{Float64}(getproperty(f, n)) for n in names)
    return Vector{Int}(f.year), M, collect(names)
end
