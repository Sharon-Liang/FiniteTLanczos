module myED

using SparseArrays
using LinearAlgebra
using QuantumLattices
using Arpack
using Random; Random.seed!()
using StatsFuns

export pauli, eye, random_init
export icgs, itFOLM

export model, TFIsing
export ave_sx
export partitian, free_energy, thermal_average, correlation2time

export FED

include("utilities.jl")
include("setups.jl")
include("methods.jl")
include("PhysicalObservables.jl")
include("exact.jl")


end # module
