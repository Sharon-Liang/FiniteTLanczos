#module PhysicalObservables

"""Thermal dynamic quantities"""
function partitian(β::Real, A::FED)
    e = A.val .- A.val[1]
    z = exp.(-β * e)|> sum
    return z
end

function partitian(β::Real, A::FTLM)
    n = size(A.vec)[1] / A.R
    z = 0.0
    for r = 1: A.R, j = 1: A.M
        e = A.val[j,r] - A.val[1,r]
        fac = (A.initv[:,r]' * A.vec[:,j,r]) * 
              (A.vec[:,j,r]' * A.initv[:,r])
        z += exp(-β * e) * fac
    end
    return z * n
end

function partitian(β::Real, A::OFTLM)
    n = (size(A.vec)[1] - A.Ne) / A.R
    z = 0.0
    for r = 1: A.R, j = 1: A.M
        e = A.val[j,r] - A.eval[1]
        fac = (A.initv[:,r]' * A.vec[:,j,r]) * 
              (A.vec[:,j,r]' * A.initv[:,r])
        z += exp(-β * e) * fac
    end
    z = z * n
    e = A.eval .- A.eval[1]
    z0 = exp.(-β * e) |> sum
    return z + z0
end

function free_energy(β::Real, A::ED)
    f = partitian(β, A) |> log
    return -f/(β * A.m.L)
end


function energy(β::Real, A::FED) 
    e = A.val .- A.val[1]
    res = e .* exp.(-β*e) |> sum
    return res/(partitian(β,A) * A.m.L)
end

function energy(β::Real, A::FTLM) 
    n = size(A.vec)[1] / A.R
    res = 0.0
    for r = 1: A.R, j = 1: A.M
        e = A.val[j,r] - A.val[1,r]
        fac = (A.initv[:,r]' * A.vec[:,j,r]) * 
              (A.vec[:,j,r]' * A.initv[:,r])
        res += e * exp(-β * e) * fac
    end
    return res * n /(partitian(β,A) * A.m.L)
end

function energy(β::Real, A::OFTLM) 
    n = (size(A.vec)[1] - A.Ne) / A.R
    res = 0.0
    for r = 1: A.R, j = 1: A.M
        e = A.val[j,r] - A.eval[1]
        fac = (A.initv[:,r]' * A.vec[:,j,r]) * 
              (A.vec[:,j,r]' * A.initv[:,r])
        res += e * exp(-β * e) * fac
    end
    res = res * n
    e = A.eval .- A.eval[1]
    res0 = e .* exp.(-β*e) |> sum
    return (res+res0)/(partitian(β,A) * A.m.L)
end

function specific_heat(β::Real, A::FED)
    e = A.val .- A.val[1]
    c = e .* e .* exp.(-β*e) |> sum
    c = c / partitian(β,A) - energy(β,A)^2
    return c * β^2 / A.m.L
end

function specific_heat(β::Real, A::FTLM) 
    n = size(A.vec)[1] / A.R
    res = 0.0
    for r = 1: A.R, j = 1: A.M
        e = A.val[j,r] - A.val[1,r]
        fac = (A.initv[:,r]' * A.vec[:,j,r]) * 
              (A.vec[:,j,r]' * A.initv[:,r])
        res += e^2 * exp(-β * e) * fac
    end
    res = res * n / partitian(β,A) - energy(β,A)^2
    return res * β^2 / A.m.L
end

function specific_heat(β::Real, A::OFTLM) 
    n = (size(A.vec)[1] - A.Ne) / A.R
    res = 0.0
    for r = 1: A.R, j = 1: A.M
        e = A.val[j,r] - A.eval[1]
        fac = (A.initv[:,r]' * A.vec[:,j,r]) * 
              (A.vec[:,j,r]' * A.initv[:,r])
        res += e^2 * exp(-β * e) * fac
    end
    res = res * n
    e = A.eval .- A.eval[1]
    res0 = e.^2 .* exp.(-β*e) |> sum
    res = (res + res0)/partitian(β,A) - energy(β,A)^2
    return res * β^2 / A.m.L
end

function entropy(β::Real, A::ED)
    return energy(β,A)*β - β * free_energy(β,A) 
end


"""Thermal average"""
function thermal_average(β::Real, O::AbstractMatrix, A::FED)
    d1 = size(O)[1]; d2 = length(A.val)
    d = d2/d1 |> Integer
    O = O ⊗ eye(d)
    O = A.vec' * O * A.vec
    e = A.val .- A.val[1]
    res = exp.(-β*e) .* diag(O) |> sum
    z = partitian(β, A)
    return res/z
end

function thermal_average(β::Real, O::AbstractMatrix, A::FTLM)
    d1 = size(O)[1]; d2 = size(A.vec)[1]
    d = d2/d1 |> Integer
    O = O ⊗ eye(d)
    n = size(A.vec)[1] / A.R
    res = 0.0
    for r = 1: A.R, j = 1:A.M
        e = A.val[j,r] - A.val[1,r]
        fac = (A.initv[:,r]' * A.vec[:,j,r]) * 
              (A.vec[:,j,r]' * O * A.initv[:,r])
        res += exp(-β*e) * fac
    end
    return res*n/partitian(β, A)
end

function thermal_average(β::Real, O::AbstractMatrix, A::OFTLM)
    d1 = size(O)[1]; d2 = size(A.vec)[1]
    d = d2/d1 |> Integer
    O = O ⊗ eye(d)
    n = (size(A.vec)[1] - A.Ne) / A.R
    res = 0.0
    for r = 1: A.R, j = 1: A.M
        e = A.val[j,r] - A.eval[1]
        fac = (A.initv[:,r]' * A.vec[:,j,r]) * 
              (A.vec[:,j,r]' * O * A.initv[:,r])
        res += exp(-β*e)*fac
    end
    res = res * n
    O = A.evec' * O * A.evec
    e = A.eval .- A.eval[1]
    res0 = exp.(-β*e) .* diag(O) |> sum
    return (res+res0)/partitian(β, A)
end

"""Correlations and Structure factor"""
function c_average(N::Integer, m::model; cf::Function, mf::Function)
    ave = 0.0
    for i = 1:N
        h = mf(m)
        ave += cf(m)
    end
    return ave/N
end

# ----------------------------------------------------------------------#
# ----------------------------------------------------------------------#
function correlation2time(τ::Real, β::Real,
    O1::T, O2::T, A::FED) where T<:AbstractMatrix
    d1 = size(O1)[1]; d2 = length(A.val)
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

function correlation2time(τ::Real, β::Real,
        O1::T, O2::T, A::FTLM) where T<:AbstractMatrix
    d1 = size(O1)[1]; d2 = size(A.vec)[1]
    d = d2/d1 |> Integer
    O1 = O1 ⊗ eye(d) ; O2 = O2 ⊗ eye(d)
    n = size(A.vec)[1] / A.R
    res = 0.0
    for r = 1: A.R
        ei = A.val[:,r] .- A.val[1,r]
        #2nd Lanczos procedure
        v0 = O2 * A.initv[:,r]
        v0 = v0 / norm(v0)
        T2, Q2 = itFOLM(A.m.ham, v0, nev = A.M)
        ej, v = eigen(T2)
        ej = ej .- ej[1]
        vec = -Q2 * v
        for i = 1: A.M, j = 1: A.M
            fac = (A.initv[:,r]' * A.vec[:,i,r]) * 
                (A.vec[:,i,r]' * O1 * vec[:,j]) *
                (vec[:,j]' * O2 * A.initv[:,r])
            res += exp(-β*ei[i]+τ*ei[i]-τ*ej[j]) * fac
        end
    end
    return res*n/partitian(β,A)
end

"""
function correlation2time(τ::Real, β::Real,
    O1::T, O2::T, A::OFTLM) where T<:AbstractMatrix
    d1 = size(O1)[1]; d2 = size(A.vec)[1]
    d = d2/d1 |> Integer
    O1 = O1 ⊗ eye(d) ; O2 = O2 ⊗ eye(d)
    n = (size(A.vec)[1] - A.Ne) / A.R
    res = 0.0
    for r = 1: A.R
        ei = A.val[:,r] .- A.eval[1]
        #2nd Lanczos procedure
        v0 = O2 * A.initv[:,r]; v0 = v0 / norm(v0)
        v0 = icgs(v0, A.evec)
        T2, Q2 = itFOLM(A.m.ham, [A.evec v0], nev = A.M-1)
        Q2 = [v0 Q2]
        T2 = Q2' * A.m.ham * Q2
        ej, v = eigen(T2)
        ej = ej .- A.eval[1]
        vec = Q2 * v
        for i = 1:A.M, j = 1: A.M
            fac = (A.initv[:,r]' * A.vec[:,i,r]) * 
                (A.vec[:,i,r]' * O1 * vec[:,j]) *
                (vec[:,j]' * O2 * A.initv[:,r])
            res += exp(-β*ei[i]+τ*ei[i]-τ*ej[j]) * fac
        end
    end
    res = res * n
    O1 = A.evec' * O1 * A.evec
    O2 = A.evec' * O2 * A.evec
    e = A.eval .- A.eval[1]
    res0 = 0.0
    for i = 1: A.Ne, j = 1: A.Ne
        res0 += exp(-β*e[i]+τ*(e[i] - e[j])) * O1[i,j] * O2[j,i]
    end
    return (res + res0)/partitian(β,A)
end
"""

# ----------------------------------------------------------------------#
# ----------------------------------------------------------------------#
function imag_susceptibility(ω::Real, β::Real,
    O1::T, O2::T, A::FED; η::Real = 0.05) where T<:AbstractMatrix
    d1 = size(O1)[1]; d2 = length(A.val)
    d = d2/d1 |> Integer
    O1 = O1 ⊗ eye(d) ; O2 = O2 ⊗ eye(d)
    O1 = A.vec' * O1 * A.vec
    O2 = A.vec' * O2 * A.vec
    e = A.val .- A.val[1]
    num = 0.0
    for i = 1: d2, j = 1: d2
        res = exp(-β*e[i]) - exp(-β*e[j])
        res = res * O1[i,j] * O2[j,i] * delta(ω+e[i]-e[j],η)
        num += res
    end
    return  π*num/partitian(β,A)
end

function imag_susceptibility(ω::Real, β::Real,
    O1::T, O2::T, A::FTLM; η::Real = 0.05) where T<:AbstractMatrix
    d1 = size(O1)[1]; d2 = size(A.vec)[1]
    d = d2/d1 |> Integer
    O1 = O1 ⊗ eye(d) ; O2 = O2 ⊗ eye(d)
    n = size(A.vec)[1] / A.R
    res = 0.0
    for r = 1: A.R
        ei = A.val[:,r] .- A.val[1,r]
        #2nd Lanczos procedure
        v0 = O2 * A.initv[:,r]
        v0 = v0 / norm(v0)
        T2, Q2 = itFOLM(A.m.ham, v0, nev = A.M)
        ej, v = eigen(T2)
        ej = ej .- A.val[1,r]
        vec = Q2 * v
        for i = 1:A.M, j = 1: A.M
            fac = (A.initv[:,r]' * A.vec[:,i,r]) * 
                (A.vec[:,i,r]' * O1 * vec[:,j]) *
                (vec[:,j]' * O2 * A.initv[:,r])
            res += (exp(-β*ei[i])-exp(-β*ej[j]))*fac*delta(ω+ei[i]-ej[j],η)
        end
    end
    res = res * n
    return  π*res/partitian(β,A)
end

# ----------------------------------------------------------------------#
# ----------------------------------------------------------------------#
function structure_factor(ω::Real, β::Real,
    O1::T, O2::T, A::FED; η::Real = 0.05) where T<:AbstractMatrix
    d1 = size(O1)[1]; d2 = length(A.val)
    d = d2/d1 |> Integer
    O1 = O1 ⊗ eye(d) ; O2 = O2 ⊗ eye(d)
    O1 = A.vec' * O1 * A.vec
    O2 = A.vec' * O2 * A.vec
    e = A.val .- A.val[1]
    num = 0.0
    for i = 1: d2, j = i+1: d2
        num += exp(-β*e[i]) * O1[i,j] * O2[j,i] * delta(ω+e[i]-e[j],η)
    end
    return  2π * num/partitian(β,A)
end

function structure_factor(ω::Real, β::Real,
    O1::T, O2::T, A::FTLM; η::Real = 0.05) where T<:AbstractMatrix
    d1 = size(O1)[1]; d2 = size(A.vec)[1]
    d = d2/d1 |> Integer
    O1 = O1 ⊗ eye(d) ; O2 = O2 ⊗ eye(d)
    n = size(A.vec)[1] / A.R
    res = 0.0
    for r = 1: A.R
        ei = A.val[:,r] .- A.val[1,r]
        #2nd Lanczos procedure
        v0 = O2 * A.initv[:,r]
        v0 = v0 / norm(v0)
        T2, Q2 = itFOLM(A.m.ham, v0, nev = A.M)
        ej, v = eigen(T2)
        ej = ej .- A.val[1,r]
        vec = Q2 * v
        for i = 1:A.M, j = 1: A.M
            fac = (A.initv[:,r]' * A.vec[:,i,r]) * 
                (A.vec[:,i,r]' * O1 * vec[:,j]) *
                (vec[:,j]' * O2 * A.initv[:,r])
            res += exp(-β*ei[i]) *fac*delta(ω+ei[i]-ej[j],η)
        end
    end
    res = res * n
    return  2π*res/partitian(β,A)
end
#end  # modulePhysicalObservables