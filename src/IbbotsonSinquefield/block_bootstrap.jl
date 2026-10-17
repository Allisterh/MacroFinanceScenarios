# Circular (moving or stationary) block bootstrap of whole historical years: each simulated
# year is one row of the sample, so the cross-sectional correlation of a year is kept, and
# blocks of consecutive years keep the autocorrelation up to the block length.

"""
    block_bootstrap(X::AbstractMatrix, n_steps, n_scen; block_size=5, random_block=false,
                    rng=Random.default_rng()) -> Array{Float64,3}

Resample the rows (years) of the `n_obs × n_vars` sample `X` in blocks of consecutive rows
into `n_scen` paths of `n_steps` years; returns `n_vars × n_steps × n_scen`. `n_steps` may
exceed `n_obs`.

Blocks are circular: each starts at a uniform row of `1:n_obs` and wraps from the last row
to the first, so every year is drawn with the same probability and the bootstrap mean
equals the sample mean in expectation. With `random_block = false` every block has
`block_size` years (the last one of a path is cut at `n_steps`); with `random_block = true`
the lengths are `1 + Geometric(1/block_size)`, mean `block_size` (the stationary bootstrap of
Politis and Romano).
"""
function block_bootstrap(X::AbstractMatrix, n_steps::Integer, n_scen::Integer; block_size::Integer = 5,
                         random_block::Bool = false, rng::AbstractRNG = Random.default_rng())
    n_obs, n_vars = size(X)
    1 <= block_size <= n_obs ||
        throw(ArgumentError("block_size must be in 1:$n_obs (the sample length), got $block_size"))
    n_steps >= 1 || throw(ArgumentError("n_steps must be ≥ 1, got $n_steps"))
    n_scen >= 1 || throw(ArgumentError("n_scen must be ≥ 1, got $n_scen"))
    len_dist = Geometric(1 / block_size)
    out = Array{Float64,3}(undef, n_vars, n_steps, n_scen)
    for s in 1:n_scen
        t = 1
        while t <= n_steps
            start = rand(rng, 1:n_obs)
            len = random_block ? 1 + rand(rng, len_dist) : block_size
            for k in 0:min(len, n_steps - t + 1)-1
                out[:, t, s] = @view X[mod1(start + k, n_obs), :]
                t += 1
            end
        end
    end
    return out
end

"""
    block_bootstrap(ta::TimeArray, n_steps, n_scen; start_year=year(last(timestamp(ta))) + 1,
                    kwargs...) -> Scenarios

Block bootstrap of the columns of `ta` (see the matrix method for `kwargs`), as `Scenarios`
named by `colnames(ta)` whose first simulated year is `start_year`.
"""
function block_bootstrap(ta::TimeArray, n_steps::Integer, n_scen::Integer;
                         start_year::Integer = year(last(timestamp(ta))) + 1, kwargs...)
    return Scenarios(block_bootstrap(values(ta), n_steps, n_scen; kwargs...), colnames(ta), start_year)
end
