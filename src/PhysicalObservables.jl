#module PhysicalObservables
function partitian(β::Number, A::FED)
    e = A.val .- A.val[1]
    z = exp.(-β * e)|> sum
    return z
end

function free_energy(β::Number, A::FED)
    f = partitian(β, A) |> log
    return -f/(A.L*β)
end

function thermal_average(β::Number, O::AbstractMatrix, A::FED)
    d1 = size(O)[1]; d2 = length(A.val)
    d = d2/d1 |> Integer
    O = O ⊗ eye(d)
    O = A.vec' * O * A.vec
    e = A.val .- A.val[1]
    res = exp.(-β*e) .* diag(O) |> sum
    z = partitian(β, A)
    return res/z
end

"""
function correlation2time(τ::Number, β::Number,
        O1::T, O2::T, A::FED) where T<:AbstractMatrix
    d1 = size(O)[1]; d2 = length(A.val)
    d = d2/d1 |> Integer
    O1 = O1 ⊗ eye(d) ; O2 = O2 ⊗ eye(d)
    O1 = A.vec' * O1 * A.vec
    O2 = A.vec' * O2 * A.vec
    e = A.val .- A.val[1]
    num = 0.0
    for i = 1: d2, j = 1: d2
        num += exp(-β*e[i]+τ*(e[i] - e[j])) * O1[i,j] * O2[j,i]
    end
    return num/partitian(β,A)
end
"""
#end  # modulePhysicalObservables
