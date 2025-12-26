module SURD

include("Utils.jl")
using .Utils
using StatsBase

export surd

"""
    surd(target::AbstractVector, agents::Tuple)

Calculate the Synergistic, Unique, and Redundant Decomposition (SURD) of information provided by multiple agents about a target variable.

# Arguments
- `target`: A vector representing the target variable.
- `agents`: A tuple of vectors representing the agent variables.
- `log`: (optional) If true, prints the SURD results to the console.
"""
function surd(
    target::AbstractVector,
    agents::Tuple;
    log::Bool=false,
)

    # Get the dimension of the target vector
    agents_dimension = length(agents)

    # probability of full vector
    p_full = fit(Histogram, (target, agents...))
    # probability of target only
    target_edges = p_full.edges[1]
    p_target = fit(Histogram, target, target_edges)
    pdf_target = p_target.weights / sum(p_target.weights)
    # probability of agents only
    p_agents = fit(Histogram, agents, p_full.edges[2:end])

    # leak entropy --------------------------------------------------
    H_full = shannon_entropy(p_full)
    H_agents = shannon_entropy(p_agents)
    H_target = shannon_entropy(p_target)
    H_leak = H_full - H_agents
    H_leak_norm = H_leak / H_target

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
    for target_index in 1:target_edges_length-1

        MI *= 0.0# reset MI for each target bin
        ΔI = 0.0

        for dim in 1:agents_dimension
            for num in 1:length(comb[dim])
                MI[dim][num] = particular_mutual_info(p_target, p_agents, p_full, target_index, comb[dim][num])
            end
        end

        # Input I_r
        sort_num = sortperm(MI[1])
        num = copy(comb[end][end])# full set of agents
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

    # Output results
    if log
        println("- SURD Results ------------------------------")
        for dim in 2:agents_dimension
            for num in 1:length(comb[dim])
                println("I_s(", comb[dim][num], ") \t$(round(I_s[dim-1][num] / I_total, digits=3))")
            end
        end
        println("- - - - - - - - - - - -")
        for dim in 2:agents_dimension
            for num in 1:length(comb[dim])
                println("I_r(", comb[dim][num], ") \t$(round(I_r[dim-1][num] / I_total, digits=3))")
            end
        end
        println("- - - - - - - - - - - -")
        for i in 1:agents_dimension
            println("I_u(", i, ") \t\t$(round(I_u[i] / I_total, digits=3))")
        end
        println("- - - - - - - - - - - -")
        println("H_leak: \t$(round(H_leak_norm, digits=3))")
        println("----------------------------------------------")
    end

    I_u_dict = Dict()
    I_r_dict = Dict()
    I_s_dict = Dict()

    for i in 1:agents_dimension
        push!(I_u_dict, "$i" => I_u[i] / I_total)
    end
    for dim in 2:agents_dimension
        for num in 1:length(comb[dim])
            push!(I_s_dict, string(join(comb[dim][num], ",")) => I_s[dim-1][num] / I_total)
            push!(I_r_dict, string(join(comb[dim][num], ",")) => I_r[dim-1][num] / I_total)
        end
    end

    return Dict("u" => I_u_dict, "s" => I_s_dict, "r" => I_r_dict, "leak" => H_leak_norm)
end

end
