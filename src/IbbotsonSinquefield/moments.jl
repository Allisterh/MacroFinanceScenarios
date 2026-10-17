# Moment matching: shift columns so that their sample means hit given targets, e.g. to
# impose a long-run equity return on historical or simulated data.

"""
    match_means(ta::TimeArray, targets::NamedTuple) -> TimeArray
    match_means(sc::Scenarios, targets::NamedTuple) -> Scenarios

Shift each variable named in `targets` by a constant so that its sample mean equals the
target, e.g. `match_means(ta, (re = 0.05, rb = 0.02))`. Variables not in `targets` are
left unchanged. For `Scenarios` the mean is taken over all years and paths. Fails if a
target names a variable missing from the data.
"""
function match_means(ta::TimeArray, targets::NamedTuple)
    names = colnames(ta)
    check_target_names(names, targets)
    X = copy(values(ta))
    for (c, μ) in pairs(targets)
        j = findfirst(==(c), names)
        X[:, j] .+= μ - mean(view(X, :, j))
    end
    return TimeArray(timestamp(ta), X, names)
end

function match_means(sc::Scenarios, targets::NamedTuple)
    check_target_names(sc.names, targets)
    data = copy(sc.data)
    for (c, μ) in pairs(targets)
        i = findfirst(==(c), sc.names)
        data[i, :, :] .+= μ - mean(view(data, i, :, :))
    end
    return Scenarios(data, sc.names, sc.start_year)
end

function check_target_names(names, targets::NamedTuple)
    unknown = setdiff(keys(targets), names)
    isempty(unknown) || throw(ArgumentError("unknown variables: $(join(unknown, ", "))"))
end
