module SURDCairoMakieExt

using CairoMakie
using LaTeXStrings

import SURD
import SURD: surd_plot

function format_label(x)
    if length(x) == 1
        return string(x[1])
    else
        return x
    end
end

function surd_plot(result::Dict;
    threshold::Float64=0.01,
    size=nothing,
    colorset=[:lightcoral, :cyan4, :darkgoldenrod1, :gray]
)

    labels = String[]
    vals = Float64[]
    colors = Symbol[]

    member_length = 0

    # u
    for (k, v) in result["u"]
        if v > threshold
            tmp = format_label.(k)
            push!(labels, L"U_{%$tmp}")
            push!(vals, v)
            push!(colors, colorset[1])
            member_length += 1
        end
    end
    # r
    for (k, v) in result["r"]
        if v > threshold
            tmp = format_label.(k)
            push!(labels, L"R_{%$tmp}")
            push!(vals, v)
            push!(colors, colorset[2])
            member_length += 1
        end
    end

    # s
    for (k, v) in result["s"]
        if v > threshold
            tmp = format_label.(k)
            push!(labels, L"S_{%$tmp}")
            push!(vals, v)
            push!(colors, colorset[3])
            member_length += 1
        end
    end

    figure = with_theme(theme_latexfonts()) do

        fig = Figure(size=((member_length+1) * 100, 300))
        if !isnothing(size)
            fig = Figure(size=size)
        end

        ax1 = Axis(fig[1, 1],
            xticks=(1:length(labels), latexstring.(labels)),
            xticklabelrotation=π / 4,
            limits=(nothing, (0, 1)),
            title=L"\Delta I_{(\cdot)\rightarrow j}/I(Q_j^+:\mathbf{Q})"
        )

        barplot!(ax1, 1:length(vals), vals, color=colors, strokecolor=:black, strokewidth=1.0)

        # leak plot
        ax2 = Axis(fig[1, 2],
            xticks=([1], [""]),
            limits=(nothing, (0, 1)),
            title="leak"
        )

        barplot!(ax2, [1], [result["leak"]], color=colorset[4], strokecolor=:black, strokewidth=1.0)

        colsize!(fig.layout, 1, Relative(member_length*1.5 / (member_length*1.5 + 1)))
        colsize!(fig.layout, 2, Relative(1 / (member_length*1.5 + 1)))

        fig
    end
end

end
