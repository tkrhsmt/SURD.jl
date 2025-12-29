```@meta
CurrentModule = SURD
```

# SURD

Documentation for [SURD](https://github.com/tkrhsmt/SURD.jl).
SURD (Synergistic-Unique-Redundant Decomposition) is a Julia package for performing SURD decomposition on multivariate time series data.

## Installation

You can install SURD.jl using Julia's package manager. Open the Julia REPL and run the following commands (press `]` to enter the package manager):

```julia
add https://github.com/tkrhsmt/SURD.jl
```


## Usage

After installing the package, you can use it in your Julia scripts or REPL as follows:

```julia
using SURD

# Example data
Q1 = zeros(1000000)
Q2 = zeros(1000000)
Q3 = zeros(1000000)

# Initialize the time series with some random values
Q1[1] = randn()
Q2[1] = randn()
Q3[1] = randn()

# Generate synthetic time series data
for i in 2:1000000
    Q1[i] = 0.3*Q1[i-1] + sin(Q2[i-1] * Q3[i-1]) + 0.001 * randn()
    Q2[i] = 0.5 * Q2[i-1] + 0.1 * randn()
    Q3[i] = Q2[i]
end

# Perform SURD decomposition
# (Q1, Q2, Q3) -> Q2^+
result = surd(Q2[2:end], (Q1[1:end-1], Q2[1:end-1], Q3[1:end-1]); log=true)

```

## Visualization

You can visualize the results of the SURD decomposition using the following code:

```julia
using CairoMakie

fig = Figure(size=(700, 200))
ax1 = Axis(fig[1, 1], xticks=(1:11, ["U1", "U2", "U3", "R12", "R13","R23", "R123", "S12", "S13", "S23", "S123"]), aspect=5, xticklabelrotation=pi / 4, height=100, width=500)

barplot!(ax1, [1:11;],
    [
        result["u"]["1"], result["u"]["2"], result["u"]["3"],
        result["r"]["1,2"], result["r"]["1,3"], result["r"]["2,3"], result["r"]["1,2,3"],
        result["s"]["1,2"], result["s"]["1,3"], result["s"]["2,3"], result["s"]["1,2,3"],
    ],
    color=[1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
)
ylims!(ax1,(0,1))

ax2 = Axis(fig[1, 2], xticks=(1:1, ["leak"]), aspect=1/5, height=100, width=80)
barplot!(ax2, [1], [result["leak"]], color = :gray)
ylims!(ax2,(0,1))

save("surd_decomposition.png", fig)
```
