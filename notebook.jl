### A Pluto.jl notebook ###
# v0.11.12

using Markdown
using InteractiveUtils

# ╔═╡ 24f312c2-f126-11ea-37a8-d3ad904331c3
begin
	using LinearAlgebra
	n = 10
end

# ╔═╡ 46c88734-f118-11ea-363d-51e5669d04b9
using Plots

# ╔═╡ 8f168358-f10c-11ea-2ee3-17046db1da2e
md"# Percolation"

# ╔═╡ 8ea56784-f11b-11ea-2db8-3d1fd57dad9d
md"""
[Percolation](https://en.wikipedia.org/wiki/Percolation_theory) describes a phenomenon where the "connectedness" of a system changes suddenly at a certain entity density.

I came across this concept while listening to the German [Coronavirus-Update podcast](https://www.youtube.com/watch?v=Yho1g3o9EZs) with Prof. Christian Drosten, where he describes how percolation can unpredictably change the spreading of SARS-CoV-2.

To have an accessable visualisation of percolation, I've put together this little notebook.
Our system is very simple: An ``m × n`` matrix with zeros and ones where *ones represent entities and zeros represent empty spots*.
Two entities are said to be *connected* if they are at adjacent positions in the matrix (von Neumann neighborhood).
Finally, we call the whole system *connected* iff there is a path of connected ones from the top left to bottom right corner.

As an example:
This system is connected:
```math
\begin{pmatrix}
	1 & 1 & 0 & 1 \\
	0 & 1 & 1 & 0 \\
	1 & 1 & 0 & 0 \\
	1 & 0 & 0 & 0 \\
	1 & 1 & 1 & 1
\end{pmatrix}
```
and this one is not:
```math
\begin{pmatrix}
	0 & 1 & 0 & 1 \\
	0 & 0 & 1 & 0 \\
	1 & 1 & 0 & 1 \\
	0 & 0 & 0 & 1
\end{pmatrix}
```

The question to answer is: **If the proportion of ones is ``p``, what is the probability that the system is connected?**
"""

# ╔═╡ 8fbd5282-f11e-11ea-0413-7378f02e1856
md"""
The following function implements checking for such connectedness given a matrix of `Bool` values and from where to where the connection should be checked (defaults to top left to bottom right).
The algorithm is basically a breadth-first-search.
"""

# ╔═╡ f50e403a-f112-11ea-3c53-df2855ecf4f7
const von_neumann_neighbors_2 =
	( CartesianIndex(0, 1)
	, CartesianIndex(0, -1)
	, CartesianIndex(1, 0)
	, CartesianIndex(-1, 0)
	)

# ╔═╡ c21785c0-f10c-11ea-1d81-13894c5629f5
function check_connection(grid::AbstractMatrix{Bool};
		start=first(CartesianIndices(grid)),
		stop=last(CartesianIndices(grid)))
	
	idcs = CartesianIndices(grid)
	first_idx, last_idx = first(idcs), last(idcs)
	
	all_reached = Set{eltype(idcs)}()
	recently_reached = Set{eltype(idcs)}()
	
	if grid[start]
		push!(all_reached, start)
		push!(recently_reached, start)
	else
		return false
	end
	
	while !isempty(recently_reached)
		new_reached = Set{eltype(idcs)}()
		for idx in recently_reached
			for d in von_neumann_neighbors_2
				new_idx = max(min(idx + d, last_idx), first_idx)
				new_idx ∈ all_reached && continue
				if grid[new_idx]
					push!(all_reached, new_idx)
					push!(new_reached, new_idx)
				end
			end
		end
		recently_reached = new_reached
	end
	
	stop ∈ all_reached
end

# ╔═╡ d4e64fb2-f126-11ea-209c-6fb41a0b3e92
md"""
We can test this function by trying out some obvious examples:
"""

# ╔═╡ 2c67936a-f127-11ea-2be0-07067fa7e8d9
zigzag = Bidiagonal(trues(n), trues(n-1), :U)

# ╔═╡ 33f3e9e4-f127-11ea-39dc-9d162ef8f6cc
diagonal = Diagonal(trues(n))

# ╔═╡ 02af7e2a-f127-11ea-3e69-c3ff02c82117
check_connection(zigzag),
check_connection(zigzag; start=CartesianIndex(1, n), stop=CartesianIndex(n, 1)),
check_connection(diagonal)

# ╔═╡ 3e3ea022-f11f-11ea-1f1a-1bcde8a9cabd
md"""
Since we want to make a probabilistic statement, we should check if there is a connection through the system for multiple systems with a given proportion ``p`` of ones.

The following function estimates the probability of having a complete connection:
"""

# ╔═╡ 93e403d0-f117-11ea-19b2-e17d1bb4926f
function conn_prob(p; n=10, iter=100)
	grid = falses(n, n)
	counter = 0
	for i in 1:iter
		for idx in eachindex(grid)
			grid[idx] = rand() < p
		end
		counter += check_connection(grid)
	end
	counter / iter
end

# ╔═╡ 7847735e-f118-11ea-0f69-8149367dbd65
begin
	ps = 0:0.05:1
	plot(ps,
		conn_prob.(ps; n=100, iter=1000),
		ylabel="prob. of connectedness",
		xlabel="\$p\$",
		label=nothing)
end

# ╔═╡ f33aea44-f11f-11ea-05ea-4d6042701014
md"""
Now, let's have a look at the result!
Going in steps of $(step(ps)) through the range from $(minimum(ps)) to $(maximum(ps)) for ``p``, we plot the probability of a full connection.
"""

# ╔═╡ ff901008-f124-11ea-094f-6de9a8f65b77
md"""
We see that there is not at all a linear relationship but the connectedness probability suddenly increases at a very specific value for ``p``.
"""

# ╔═╡ 61b4d778-f125-11ea-153f-75fa405de8d1
md"""
Obviously, the existence of a continuous path between the upper left and lower right corner is just a proxy for how connected the system is in general.
Try using other values for `start` and `stop` in the call of `check_connection` in `conn_prob` to see if this choice matters.
"""

# ╔═╡ Cell order:
# ╟─8f168358-f10c-11ea-2ee3-17046db1da2e
# ╟─8ea56784-f11b-11ea-2db8-3d1fd57dad9d
# ╟─8fbd5282-f11e-11ea-0413-7378f02e1856
# ╠═c21785c0-f10c-11ea-1d81-13894c5629f5
# ╠═f50e403a-f112-11ea-3c53-df2855ecf4f7
# ╟─d4e64fb2-f126-11ea-209c-6fb41a0b3e92
# ╠═24f312c2-f126-11ea-37a8-d3ad904331c3
# ╠═2c67936a-f127-11ea-2be0-07067fa7e8d9
# ╠═33f3e9e4-f127-11ea-39dc-9d162ef8f6cc
# ╠═02af7e2a-f127-11ea-3e69-c3ff02c82117
# ╟─3e3ea022-f11f-11ea-1f1a-1bcde8a9cabd
# ╠═93e403d0-f117-11ea-19b2-e17d1bb4926f
# ╟─f33aea44-f11f-11ea-05ea-4d6042701014
# ╠═46c88734-f118-11ea-363d-51e5669d04b9
# ╠═7847735e-f118-11ea-0f69-8149367dbd65
# ╟─ff901008-f124-11ea-094f-6de9a8f65b77
# ╟─61b4d778-f125-11ea-153f-75fa405de8d1
