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
estimation windows, the starting point `y0` (state `[π, r, tp, erp]` in `last_year`) and
`model_params`, one entry per state variable. π, r and tp hold their AR(1) fits
(`n, c, se_c, phi, se_phi, sigma, r2`); erp holds only `c = erp_target` and `phi = 0`. Each
value is stored once: `load_issm_params` builds the VAR intercepts and lags from the `c` and
`phi` entries. The AR window is read from the fits, which must share it.
"""
function save_issm_params(path::AbstractString; ar_π::AR1, ar_r::AR1, ar_tp::AR1,
                          erp_target::Real, boot_period::Tuple{Integer,Integer},
                          y0::AbstractVector{<:Real}, last_year::Integer)
    fits = (π = ar_π, r = ar_r, tp = ar_tp)
    allequal(f.period for f in fits) ||
        throw(ArgumentError("the AR(1) fits have different periods: $([f.period for f in fits])"))
    length(y0) == length(ISSM_NAMES) ||
        throw(DimensionMismatch("y0 has length $(length(y0)), expected $(length(ISSM_NAMES))"))
    yaml_map(ks, vs) = "{" * join(("$k: $(yaml_value(v))" for (k, v) in zip(ks, vs)), ", ") * "}"

    open(path, "w") do io
        println(io, "# Ibbotson–Sinquefield specification, decimal log units (save_issm_params)")
        println(io, "ar_period: ", yaml_list(ar_π.period))
        println(io, "boot_period: ", yaml_list(boot_period))
        println(io, "last_year: ", yaml_value(last_year), "  # simulation starts in last_year + 1")
        println(io, "y0: ", yaml_map(ISSM_NAMES, y0), "  # state in last_year")
        println(io, "model_params:  # y_t = c + phi y_{t-1} + u_t; AR(1) fits on ar_period, erp i.i.d. around c")
        for (name, f) in pairs(fits)
            println(io, "  ", name, ": ", yaml_map(AR1_FIELDS, (length(f.resid), f.c, f.se[1], f.φ, f.se[2], f.σ, f.r2)))
        end
        println(io, "  erp: ", yaml_map((:c, :phi), (erp_target, 0.0)))
    end
    return path
end

"""
    load_issm_params(path) -> NamedTuple

Read a specification written by `save_issm_params`. Fields: `names` (state variables),
`c`, `A` (diagonal lag matrix), `y0`, `start_year = last_year + 1`, `last_year`,
`ar_period`, `boot_period`, `erp_target` (the erp intercept) and `model_params`, a vector of
NamedTuples `(variable, n, c, se_c, phi, se_phi, sigma, r2)`, one per state variable (a
Tables.jl table); the AR(1) statistics of erp are `missing`.
"""
function load_issm_params(path::AbstractString)
    d = YAML.load_file(path)
    names = copy(ISSM_NAMES)
    model_params = [(; variable = String(v), (k => get(d["model_params"][String(v)], String(k), missing)
                                               for k in AR1_FIELDS)...)
                    for v in names]
    c = Float64[m.c for m in model_params]
    A = Matrix{Float64}(Diagonal(Float64[m.phi for m in model_params]))
    y0 = [Float64(d["y0"][String(n)]) for n in names]
    period(x) = (Int(x[1]), Int(x[2]))
    return (; names, c, A, y0, start_year = Int(d["last_year"]) + 1,
            last_year = Int(d["last_year"]), ar_period = period(d["ar_period"]),
            boot_period = period(d["boot_period"]), erp_target = c[end], model_params)
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
