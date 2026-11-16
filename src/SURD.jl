module SURD

include("Utils.jl")
using .Utils
import StatsBase
using StatsBase: fit, mean
using NearestNeighbors: Chebyshev, KDTree, inrange!, knn, nn
using SpecialFunctions: digamma

export surd
export surd_plot

function surd_plot end

abstract type method end

"""
    SURD.Histogram()
    SURD.Histogram(nbins)
    SURD.Histogram(; nbins=nothing)

Estimate the decomposition using a joint histogram. By default, StatsBase
chooses the bins automatically. Pass a positive integer `nbins` to use exactly
that many bins along each dimension, including the target and every agent.
"""
Base.@kwdef struct Histogram <: method
    nbins::Union{Nothing,Int} = nothing

    function Histogram(nbins::Union{Nothing,Integer})
        isnothing(nbins) || nbins > 0 || throw(ArgumentError("nbins must be positive"))
        return new(nbins)
    end
end

function histogram_edges(series, nbins)
    lo, hi = extrema(series)
    isfinite(lo) && isfinite(hi) || throw(ArgumentError("histogram observations must be finite"))
    if lo == hi
        # Give a constant series a nonzero span.
        lo, hi = min(lo, 0), max(hi, 1)
    end
    lower, upper = promote(float(lo), float(hi))
    edges = collect(range(lower, upper; length=nbins + 1))
    if all(i -> edges[i] < edges[i + 1], 1:nbins)
        # Round outwards to include the data in the left-closed bins.
        edges[1] = lower > lo ? prevfloat(lower) : lower
        edges[end] = nextfloat(upper)
        return edges
    end

    # Exact arithmetic keeps the edges distinct when floating-point edges coincide.
    exact_lo, exact_hi = Rational{BigInt}.((lo, hi))
    width = (exact_hi - exact_lo) / nbins
    exact_edges = collect(range(exact_lo; step=width, length=nbins + 1))
    exact_edges[end] += width / 2
    return exact_edges
end

"""
    surd(target, agents, SURD.Histogram(; nbins=nothing); log=false)

Decompose the information provided by the tuple of scalar series `agents`
about the scalar series `target`. Series must be aligned and have equal length.

Returns dictionaries `"u"`, `"r"`, and `"s"` of component fractions. If there
is no information to apportion, all fractions are zero. The `"leak"` is
the conditional entropy divided by the target entropy, or `missing` if
the target entropy is zero.
"""
function surd(
    target::AbstractVector,
    agents::Tuple,
    method::Histogram;
    log::Bool=false,
)

    isempty(target) && throw(ArgumentError("target must contain at least one observation"))
    isempty(agents) && throw(ArgumentError("at least one agent is required"))
    all(a -> a isa AbstractVector && length(a) == length(target), agents) ||
        throw(ArgumentError("histogram agents must be vectors with the same length as target"))

    # Get the dimension of the target vector
    agents_dimension = length(agents)

    # probability of full vector
    data = (target, agents...)
    nbins = method.nbins
    p_full = if isnothing(nbins)
        fit(StatsBase.Histogram, data)
    else
        edges = map(series -> histogram_edges(series, nbins), data)
        fit(StatsBase.Histogram, data, edges)
    end
    # probability of target only
    target_edges = p_full.edges[1]
    p_target = fit(StatsBase.Histogram, target, target_edges)
    pdf_target = p_target.weights / sum(p_target.weights)
    # probability of agents only
    p_agents = fit(StatsBase.Histogram, agents, p_full.edges[2:end])

    # leak entropy --------------------------------------------------
    H_full = shannon_entropy(p_full)
    H_agents = shannon_entropy(p_agents)
    H_target = shannon_entropy(p_target)
    H_leak = H_full - H_agents
    H_leak_norm = H_target > 0 ? H_leak / H_target : missing

    # SURD ----------------------------------------------------------
    # Generate all possible combinations of agents
    comb = prob_combination(agents_dimension)

    # Initialize SURD value
    MI = map(v -> map(u -> 0.0, v), comb)
    I_s = map(v -> map(u -> 0.0, v), comb[2:end])
    I_u = zeros(agents_dimension)
    I_r = map(v -> map(u -> 0.0, v), comb[2:end])
    I_total = 0.0

    target_edges_length = length(target_edges)
    for target_index in 1:(target_edges_length-1)

        iszero(pdf_target[target_index]) && continue
        ΔI = 0.0

        for dim in 1:agents_dimension
            for num in 1:length(comb[dim])
                MI[dim][num] = max(0.0, particular_mutual_info(
                    p_target, p_agents, p_full, target_index, comb[dim][num],
                ))
            end
        end

        # Input I_r
        sort_num = sortperm(MI[1])
        num = copy(comb[end][end]) # full set of agents
        idx = 1
        for i in agents_dimension:-1:2
            num2 = findfirst(==(num), comb[i])
            I_r[i-1][num2] += (MI[1][sort_num[idx]] - ΔI) * pdf_target[target_index]
            ΔI = MI[1][sort_num[idx]]

            num = setdiff(num, sort_num[idx])
            idx += 1
        end

        # Input I_u
        I_u[num[1]] += (MI[1][sort_num[idx]] - ΔI) * pdf_target[target_index]
        ΔI = MI[1][sort_num[idx]]

        # Input I_s
        for dim in 2:agents_dimension
            sort_num = sortperm(MI[dim])

            for num in 1:length(comb[dim])
                num2 = sort_num[num]
                tmp = (MI[dim][num2] - ΔI) * pdf_target[target_index]
                if tmp > 0.0
                    I_s[dim-1][num2] += tmp
                    ΔI = MI[dim][num2]
                end
            end
        end

        I_total += ΔI * pdf_target[target_index]
    end

    normalizer = I_total > 0 ? I_total : 1.0

    # Output results
    if log
        println("- SURD Results ------------------------------")
        for dim in 2:agents_dimension
            for num in 1:length(comb[dim])
                println("I_s(", comb[dim][num], ") \t$(round(I_s[dim-1][num] / normalizer, digits=3))")
            end
        end
        println("- - - - - - - - - - - -")
        for dim in 2:agents_dimension
            for num in 1:length(comb[dim])
                println("I_r(", comb[dim][num], ") \t$(round(I_r[dim-1][num] / normalizer, digits=3))")
            end
        end
        println("- - - - - - - - - - - -")
        for i in 1:agents_dimension
            println("I_u(", i, ") \t\t$(round(I_u[i] / normalizer, digits=3))")
        end
        println("- - - - - - - - - - - -")
        println("H_leak: \t", ismissing(H_leak_norm) ? missing : round(H_leak_norm, digits=3))
        println("----------------------------------------------")
    end

    I_u_dict = Dict{String,Float64}()
    I_r_dict = Dict{String,Float64}()
    I_s_dict = Dict{String,Float64}()

    for i in 1:agents_dimension
        push!(I_u_dict, "$i" => I_u[i] / normalizer)
    end
    for dim in 2:agents_dimension
        for num in 1:length(comb[dim])
            push!(I_s_dict, string(join(comb[dim][num], ",")) => I_s[dim-1][num] / normalizer)
            push!(I_r_dict, string(join(comb[dim][num], ",")) => I_r[dim-1][num] / normalizer)
        end
    end

    return Dict("u" => I_u_dict, "s" => I_s_dict, "r" => I_r_dict, "leak" => H_leak_norm)
end

"""
    SURD.kNN(k=5)
    SURD.kNN(; k=5)

Estimate the continuous decomposition using a shared radius from the kth nearest
other observation. Requires `1 <= k < N`, where `N` is the number of observations.
"""
Base.@kwdef struct kNN <: method
    k::Int = 5
end

"""
    surd(target, agents, SURD.kNN(k); log=false)

Decompose information for aligned real-valued observations. Each series is a
vector for scalar observations or a matrix with one observation per column
for vector-valued observations. The tuple `agents` has one series per agent.

Returns component fractions in dictionaries `"u"`, `"r"`, and `"s"`, and
the remaining entropy power fraction `"leak"`. The component fractions
sum to one when their estimated total is positive. The leak is
`exp(-2 * I_total / d_Y)`, where `I_total` is the total before normalization
in nats and `d_Y` is the target dimension. When the estimated total is zero,
all component fractions are zero and the leak is one.

The neighborhood in joint space excludes the queried observation. Its radius
is reused for strict counts in projected spaces and for averaging local
scores over other target observations. A zero joint radius raises an
`ArgumentError` because the continuous density estimate needs positive volume.
"""
function surd(
    target::AbstractVecOrMat,
    agents::Tuple,
    method::kNN;
    log::Bool=false,
)

    agents_dimension = length(agents)
    agents_dimension > 0 || throw(ArgumentError("at least one agent is required"))
    points_target = point_matrix(target)
    points_agents = map(point_matrix, agents)
    N = size(points_target, 2)
    k = method.k
    1 <= k < N || throw(ArgumentError("k must satisfy 1 <= k < number of observations"))
    all(points -> size(points, 2) == N, points_agents) ||
        throw(ArgumentError("target and all agents must have the same number of observations"))

    points_agents_target_full = vcat(points_agents..., points_target)
    kdtree_target = KDTree(points_target, Chebyshev())
    kdtree_agents_target_full = KDTree(points_agents_target_full, Chebyshev())

    # Step 1: Find the joint radius and the Y neighbors for each observation.
    ϵ = zeros(N)
    neighbors_target = Vector{Vector{Int}}(undef, N)
    n_Y = zeros(Int, N)
    ψ_N = digamma(N)
    # Ignore roundoff from the four digamma terms before normalizing.
    zero_tol = 4 * eps(ψ_N)

    for i in 1:N
        point_target = view(points_target, :, i)
        point_agents_target = view(points_agents_target_full, :, i)
        skip_self = j -> j == i
        _, distances = knn(kdtree_agents_target_full, point_agents_target, k, false, skip_self)
        ϵ[i] = maximum(distances)
        isfinite(ϵ[i]) && ϵ[i] > 0 ||
            throw(ArgumentError("joint neighborhood radius for observation $i must be positive and finite"))

        idxs = inrange!(Int[], kdtree_target, point_target, prevfloat(ϵ[i]), false, skip_self)
        n_Y[i] = length(idxs)

        if isempty(idxs)
            _, nearest_distance = nn(kdtree_target, point_target, skip_self)
            tied = inrange!(Int[], kdtree_target, point_target, nearest_distance, false, skip_self)
            push!(idxs, minimum(tied))
        end
        neighbors_target[i] = idxs
    end

    # Generate all possible combinations of agents.
    comb = prob_combination(agents_dimension)
    comb_offsets = cumsum([0; length.(comb)])
    ĩ = zeros(comb_offsets[end], N)

    # Count points in each projection with the same joint radius.
    for dim in 1:agents_dimension
        for num in 1:length(comb[dim])
            subset = comb[dim][num]
            I = comb_offsets[dim] + num
            points_agents_comb = vcat(points_agents[subset]...)
            kdtree_agents = KDTree(points_agents_comb, Chebyshev())
            points_agents_target_comb = dim == agents_dimension ?
                points_agents_target_full : vcat(points_agents_comb, points_target)
            kdtree_agents_target = dim == agents_dimension ?
                kdtree_agents_target_full : KDTree(points_agents_target_comb, Chebyshev())

            for i in 1:N
                point_agents = view(points_agents_comb, :, i)
                point_agents_target = view(points_agents_target_comb, :, i)
                n_A = strict_neighborcount(kdtree_agents, point_agents, ϵ[i], i)
                n_YA = strict_neighborcount(kdtree_agents_target, point_agents_target, ϵ[i], i)

                ĩ[I, i] = ψ_N - digamma(n_Y[i] + 1) - digamma(n_A + 1) + digamma(n_YA + 1)
            end
        end
    end

    # Initialize SURD values.
    MI = map(v -> map(u -> 0.0, v), comb)
    I_s = map(v -> map(u -> 0.0, v), comb[2:end])
    I_u = zeros(agents_dimension)
    I_r = map(v -> map(u -> 0.0, v), comb[2:end])
    I_total = 0.0

    # Indices of each combination's subsets from the preceding order.
    parents = [
        [findall(x -> issubset(x, subset), comb[dim-1]) for subset in comb[dim]]
        for dim in 2:agents_dimension
    ]

    for i in 1:N
        # Step 2: Average the local scores of the selected Y neighbors.
        idxs = neighbors_target[i]
        for dim in 1:agents_dimension
            for num in 1:length(comb[dim])
                I = comb_offsets[dim] + num
                MI[dim][num] = mean(view(ĩ, I, idxs))
            end
        end

        # Step 3: Enforce nonnegativity and monotonicity under inclusion.
        Tᵢ = MI[end][end] > zero_tol ? MI[end][end] : 0.0
        for dim in 1:agents_dimension
            for num in 1:length(comb[dim])
                corrected = max(0.0, MI[dim][num])
                if dim > 1
                    for parent in parents[dim-1][num]
                        corrected = max(corrected, MI[dim-1][parent])
                    end
                end
                MI[dim][num] = min(Tᵢ, corrected)
            end
        end

        # Step 4: Input I_r.
        ΔI = 0.0
        sort_num = sortperm(MI[1])
        num = copy(comb[end][end]) # full set of agents
        idx = 1
        for dim in agents_dimension:-1:2
            num2 = findfirst(==(num), comb[dim])
            I_r[dim-1][num2] += MI[1][sort_num[idx]] - ΔI
            ΔI = MI[1][sort_num[idx]]

            num = setdiff(num, sort_num[idx])
            idx += 1
        end

        # Input I_u.
        I_u[num[1]] += MI[1][sort_num[idx]] - ΔI
        ΔI = MI[1][sort_num[idx]]

        # Input I_s.
        for dim in 2:agents_dimension
            sort_num = sortperm(MI[dim])

            for num in 1:length(comb[dim])
                num2 = sort_num[num]
                tmp = MI[dim][num2] - ΔI
                if tmp > 0.0
                    I_s[dim-1][num2] += tmp
                    ΔI = MI[dim][num2]
                end
            end
        end

        I_total += Tᵢ
    end

    # Step 5: Average each accumulated component over all observations.
    I_s = map(v -> map(u -> u / N, v), I_s)
    I_u = map(u -> u / N, I_u)
    I_r = map(v -> map(u -> u / N, v), I_r)
    I_total /= N
    H_leak = exp(-2.0 * I_total / size(points_target, 1))
    normalizer = I_total > 0 ? I_total : 1.0

    # Output results
    if log
        println("- SURD Results ------------------------------")
        for dim in 2:agents_dimension
            for num in 1:length(comb[dim])
                println("I_s(", comb[dim][num], ") \t$(round(I_s[dim-1][num] / normalizer, digits=3))")
            end
        end
        println("- - - - - - - - - - - -")
        for dim in 2:agents_dimension
            for num in 1:length(comb[dim])
                println("I_r(", comb[dim][num], ") \t$(round(I_r[dim-1][num] / normalizer, digits=3))")
            end
        end
        println("- - - - - - - - - - - -")
        for i in 1:agents_dimension
            println("I_u(", i, ") \t\t$(round(I_u[i] / normalizer, digits=3))")
        end
        println("- - - - - - - - - - - -")
        println("H_leak: \t", round(H_leak, digits=3))
        println("----------------------------------------------")
    end

    I_u_dict = Dict{String,Float64}()
    I_r_dict = Dict{String,Float64}()
    I_s_dict = Dict{String,Float64}()

    for i in 1:agents_dimension
        push!(I_u_dict, "$i" => I_u[i] / normalizer)
    end
    for dim in 2:agents_dimension
        for num in 1:length(comb[dim])
            push!(I_s_dict, string(join(comb[dim][num], ",")) => I_s[dim-1][num] / normalizer)
            push!(I_r_dict, string(join(comb[dim][num], ",")) => I_r[dim-1][num] / normalizer)
        end
    end

    return Dict("u" => I_u_dict, "s" => I_s_dict, "r" => I_r_dict, "leak" => H_leak)
end

function point_matrix(series)
    points = if series isa AbstractVector{<:Real}
        reshape(collect(Float64, series), 1, :)
    elseif series isa AbstractMatrix{<:Real}
        Matrix{Float64}(series)
    else
        throw(ArgumentError("each series must be a real vector or a real matrix with observations in columns"))
    end
    size(points, 1) > 0 || throw(ArgumentError("each series must have at least one coordinate"))
    all(isfinite, points) || throw(ArgumentError("observations must be finite"))
    return points
end

end
