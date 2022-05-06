#module utilities
import Base: kron

"""
Pauli matrices
"""
function pauli(symbol::Symbol)
    if symbol==:x return [0. 1.; 1. 0.] |> sparse
    elseif symbol==:y return [0. -1im; 1im 0.] |> sparse
    elseif symbol==:iy return [0. 1.; -1. 0.] |> sparse
    elseif symbol==:z return [1. 0.; 0. -1.] |> sparse
    elseif symbol==:+ return [0. 1.; 0. 0.] |> sparse
    elseif symbol==:- return [0. 0.; 1. 0.] |> sparse
    else
        error("The input should be :x,:y,:z,:+,:-, :iy.")
    end
end

"""
identity matrix
"""
function eye(dim::Integer)
    Matrix{Float64}(I,dim,dim)
end

function eye(T::DataType, dim::Integer)
    Matrix{T}(I,dim,dim)
end

function ⊗(A::AbstractArray, B::AbstractArray)
    kron(A,B)
end

"""
Gaussian approximate delta function
"""
function delta(x::Real, η::Real)
    num = η
    den = (x^2 + η^2) * π
    return num/den
end





#end #module utilities
