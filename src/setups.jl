#module setups
struct model
    ham::AbstractMatrix
    L::Integer #number of sites
end

function TFIsing(J::Real, Γ::Real; L::Integer=10, PBC::Bool=true)
    x = pauli('x'); z = pauli('z')
    hloc, hint = spzeros(2^L, 2^L), spzeros(2^L, 2^L)
    for i = 1:L
        hloc += -Γ * eye(2^(i-1)) ⊗ x ⊗ eye(2^(L-i))
        j = i + 1
        if j <= L
            hint += -J* eye(2^(i-1)) ⊗ z ⊗ z ⊗ eye(2^(L-i-1))
        else
            PBC ? hint += -J * z ⊗ eye(2^(L-2)) ⊗ z :
                  hint += spzeros(2^L, 2^L)
        end
    end
    return model(hloc + hint, L)
end
#end  # module setups
