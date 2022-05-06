#module utilities
import Base: kron

function pauli(char::Char)
    if char=='x' return [0. 1.; 1. 0.] |> sparse
    elseif char=='y' return [0. -1im; 1im 0.] |> sparse
    elseif char=='z' return [1. 0.; 0. -1.] |> sparse
    elseif char=='+' return [0. 1.; 0. 0.] |> sparse
    elseif char=='-' return [0. 0.; 1. 0.] |> sparse
    else
        @error "The input should be 'x','y','z','+','-'."
    end
end

function eye(N::Integer)
    Matrix(1.0I, N, N) |> sparse
end

function ⊗(A::AbstractArray, B::AbstractArray)
    kron(A,B)
end

function delta(x::Real, η::Real)
    num = η
    den = (x^2 + η^2) * π
    return num/den
end





#end  # module utilities
