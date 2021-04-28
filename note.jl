### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ 12d33584-a7df-11eb-1a94-919a257a38d3
begin
	using Pkg
	Pkg.activate("./")
end

# ╔═╡ bedb18fc-0e38-478b-85d0-a7dd7b6bb73f
using FiniteTLanczos

# ╔═╡ 888382d5-df98-4b9e-880b-d1478b14b019
begin
	using Printf
	using Plots
end

# ╔═╡ e4b79d90-80c8-4d89-9b04-fd3d606cd540
len = 3

# ╔═╡ cdd42d5a-2756-4afa-9fa8-65e2b656a3aa
β = 20

# ╔═╡ 2ac30069-a9a9-4825-b987-a62049fe29f6
sample = 100

# ╔═╡ bfbd3f40-4845-4b09-85fe-7b9b6d554677
h = TFIsing(1.,1.,L=len)

# ╔═╡ 1000f0e9-9061-4d91-9864-2b2a3ef3723f
C = OFTLM(h, R = sample)

# ╔═╡ 7472ecee-fd51-4dd6-bd4f-9ae3d4cbb201
B = FTLM(h, R = sample)

# ╔═╡ 1beef2da-7734-4bbf-9e9e-948cf57035ac
fidelity(B.m.L, sample)

# ╔═╡ 5035ddae-6a99-4bac-b192-e55c0aaa1cf4
A = FED(h)

# ╔═╡ f8d6d511-f683-4b55-b04d-8e6f883f2cc3
begin
	z = pauli('z')
	x = pauli('x')
	i2 = eye(2)
end

# ╔═╡ 6bbfc44c-b587-4049-94ca-8571719beaf9
correlation2time(0,β,z,z,B)

# ╔═╡ d0d138cb-e605-4585-8ac3-ea229d84ddb0
correlation2time(β,β,z,z,B)

# ╔═╡ Cell order:
# ╠═bedb18fc-0e38-478b-85d0-a7dd7b6bb73f
# ╠═e4b79d90-80c8-4d89-9b04-fd3d606cd540
# ╠═cdd42d5a-2756-4afa-9fa8-65e2b656a3aa
# ╠═2ac30069-a9a9-4825-b987-a62049fe29f6
# ╠═1beef2da-7734-4bbf-9e9e-948cf57035ac
# ╠═6bbfc44c-b587-4049-94ca-8571719beaf9
# ╠═d0d138cb-e605-4585-8ac3-ea229d84ddb0
# ╠═1000f0e9-9061-4d91-9864-2b2a3ef3723f
# ╠═7472ecee-fd51-4dd6-bd4f-9ae3d4cbb201
# ╠═5035ddae-6a99-4bac-b192-e55c0aaa1cf4
# ╠═bfbd3f40-4845-4b09-85fe-7b9b6d554677
# ╠═f8d6d511-f683-4b55-b04d-8e6f883f2cc3
# ╠═888382d5-df98-4b9e-880b-d1478b14b019
# ╠═12d33584-a7df-11eb-1a94-919a257a38d3
