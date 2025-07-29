# GraphsOfConvexSets.jl

Julia implementation of search algorithms for graphs of convex sets.
Based on the paper: "[Shortest Paths in Graphs of Convex Sets](https://arxiv.org/abs/2101.11565)" by Tobia Marcucci, Jack Umenberger, Pablo A. Parrilo, and Russ Tedrake.

## Installation

Run the following lines in Julia:
```julia
import Pkg
Pkg.add("GraphsOfConvexSets")
```
This package requires optimization solvers to be installed.
A free option is Pajarito, which can be installed as
```julia
Pkg.add("Pajarito")
```

## Usage example

See [the following shortest path example](examples/shortest_path.jl)
