# GraphsOfConvexSets.jl

Julia implementation of search algorithms for graphs of convex sets.
Based on the paper: "[Shortest Paths in Graphs of Convex Sets](https://arxiv.org/abs/2101.11565)" by Tobia Marcucci, Jack Umenberger, Pablo A. Parrilo, and Russ Tedrake.

## Installation

Run the following lines in Julia:
```
import Pkg
Pkg.add("GraphsOfConvexSets")
```
This package requires optimization solvers to be installed.
A free option is Pajarito, which can be installed as
```
Pkg.add("Pajarito")
```

## Usage example

In the following example we construct a graph of convex sets with 3 vertices and 3 edges, and we seek a shortest path from the source vertex to the tartget vertex.

```
using GraphsOfConvexSets
using Pajarito

# initialize graph
gcs = graph_of_convex_sets()

# vertices
@vertex(gcs, s) # source
@vertex(gcs, t) # target
@vertex(gcs, v) # intermediate vertex

# variables attached to the vertices
@variable(s, xs[1:2])
@variable(t, xt[1:2])
@variable(v, xv[1:2])

# convex sets attached to the vertices
@constraint(s, xs == [0, 0])
@constraint(t, xt == [1, 1])
@constraint(v, [0.9, 0] <= xv <= [1, 0.1])

# edges
@edge(gcs, st, s, t)
@edge(gcs, sv, s, v)
@edge(gcs, vt, v, t)

# convex edge-length functions
@length(st, (xt[1] - xs[1])^2 + (xt[2] - xs[2])^2)
@length(sv, (xv[1] - xs[1])^2 + (xv[2] - xs[2])^2)
@length(vt, (xt[1] - xv[1])^2 + (xt[2] - xv[2])^2)

# convex constraints on the edges
@constraint(sv, xs[2] == xv[2])

# solve shortest-path problem
path = shortest_path!(gcs, s, t, Pajarito.Optimizer)
print(path)
print(xv)
```
