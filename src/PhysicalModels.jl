#module PhysicalModels

"""
h: Hamiltonian matrix
"""
struct MODEL{Tv<:Number, Ti<:Integer}
    h::SparseMatrixCSC{Tv,Ti}
    L::Int64 #number of sites
end


"""
NN Transvers field Ising model
    H = ∑ J Zi Zj + ∑ Γ Xi
"""
function TFIsing(J::Real, Γ::Real; L::Integer=10, PBC::Bool=true)
    x = pauli(:x); z = pauli(:z)
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
    return MODEL(hloc + hint, L)
end


"""
Heisenberg model
    H = J ∑ Xi Xj + Yi Yj + Zi Zj
AFM case: J = 1.0, after unitary transformation : U=exp(iπSy) on odd sites:
    H = -0.5J ∑ (S+ S+ + S-S-) - J∑ Zi Zj
"""
function HeisenbergModel(J::Real; L::Integer=10, PBC::Bool=true)
    sp = pauli(:+); sm = pauli(:-); sz = 0.5 * pauli(:z)
    h = -J/2 * (sp ⊗ sp + sm ⊗ sm) - J * sz ⊗ sz
    hint = spzeros(2^L, 2^L)
    for i = 1:L
        j = i + 1
        if j <= L
            hint += eye(2^(i-1)) ⊗ h ⊗ eye(2^(L-i-1)) 
        else
            PBC ? hint += -J/2 * (sp⊗eye(2^(L-2))⊗sp + sm⊗eye(2^(L-2))⊗sm) - J* sz⊗eye(2^(L-2))⊗sz :
                  hint += spzeros(2^L, 2^L)
        end
    end
    return MODEL(hint, L)
end




#end #module PhysicalModels
