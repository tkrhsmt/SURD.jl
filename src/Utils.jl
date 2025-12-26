module Utils

using Combinatorics

export prob_combination
export shannon_entropy, particular_mutual_info

"""
    prob_combination(n::Int)

Generate all possible combinations of integers from 1 to n.

# Arguments
- `n::Int`: The upper limit of the range (inclusive).
"""
function prob_combination(
    n::Int,
)
    # result combinations array
    pb = []
    # 1:n range
    nrange = [1:n;]

    for i in 1:n
        pbi = []
        # all combinations of size i
        for comb in combinations(nrange, i)
            push!(pbi, comb)
        end
        push!(pb, pbi)
    end

    return pb
end


"""
shannon_entropy(p::AbstractVector)

Calculate the Shannon entropy of a probability distribution represented by histogram `p`.

# Arguments
- `p`: Histogram object representing the probability distribution.
"""
function shannon_entropy(p)

    # probability dimension
    pdim = length(p.edges)
    # data number
    ndata = sum(p.weights)
    # Shannon Entropy
    H = 0.0

    # Calculate the entropy from the histogram.
    for idx in CartesianIndices(p.weights)
        # Calculate the probability density in the bin.
        pdf = p.weights[idx] / ndata
        if pdf >= 1.0e-10 # Avoid taking the log of zero.
            # Add the contribution of the bin to the entropy.
            H -= pdf * log(pdf)
        end
    end

    return H
end

"""
particular_mutual_info(p_target, p_agents, p_full, target_index, agent_indices)

Calculate the particular mutual information for a specific target bin and a set of agent indices.

# Arguments
- `p_target`: Histogram object for the target variable.
- `p_agents`: Histogram object for the agent variables.
- `p_full`: Histogram object for the full joint distribution of target and agents.
- `target_index`: Index of the target bin for which to calculate the mutual information.
- `agent_indices`: Indices of the agents to consider in the calculation.
"""
function particular_mutual_info(p_target, p_agents, p_full, target_index, agent_indices)

    # particular mutual information
    I = 0.0
    # data number
    ndata = sum(p_full.weights)
    # agent number
    full_nagents = length(p_agents.edges)
    nagents = length(agent_indices)

    # true agents probability
    non_agent_indices = setdiff(1:full_nagents, agent_indices)
    p_agents_true = sum(p_agents.weights, dims=non_agent_indices)
    p_full_true = sum(p_full.weights, dims=non_agent_indices .+ 1)

    for idx in CartesianIndices(p_agents_true)
        # Calculate the probability density in the bin.
        px = p_target.weights[target_index] / ndata
        py = p_agents_true[idx] / ndata
        pxy = p_full_true[target_index, idx.I...] / ndata

        if pxy >= 1.0e-10 # Avoid taking the log of zero.
            # Add the contribution of the bin to the mutual information.
            I += pxy / px * log(pxy / (px * py))
        end

    end

    return I
end

end
