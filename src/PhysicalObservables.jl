#module PhysicalObservables
"""
Masubara frequencies ωn
"""
function Masubara_freq(n::Int64, β::Real; type::Symbol=:b)
    if type == :b  N = 2n
    elseif type == :f  N = 2n + 1
    else @error "type should be :b for bosons and :f for fermions" 
    end
    return N*π/β
end

"""
partitian function
"""
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


"""
free energy
"""
function free_energy(β::Real, A::LANCZOS)
    f = partitian(β, A) |> log
    return -f/(β * A.model.L)
end


"""
energy density
"""
function energy(β::Real, A::FED) 
    e = A.val .- A.val[1]
    res = e .* exp.(-β*e) |> sum
    return res/(partitian(β,A) * A.model.L)
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
    return res * n /(partitian(β,A) * A.model.L)
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
    return (res+res0)/(partitian(β,A) * A.model.L)
end


"""
specfic heat
"""
function specific_heat(β::Real, A::FED)
    e = A.val .- A.val[1]
    c = e .* e .* exp.(-β*e) |> sum
    c = c / partitian(β,A) - energy(β,A)^2
    return c * β^2 / A.model.L
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
    return res * β^2 / A.model.L
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
    return res * β^2 / A.model.L
end


"""
entropy
"""
function entropy(β::Real, A::LANCZOS)
    return energy(β,A)*β - β * free_energy(β,A) 
end


"""
Thermal average
"""
function thermal_average(β::Real, Op::AbstractMatrix, A::FED)
    d1 = size(Op)[1]; d2 = length(A.val)
    d = d2/d1 |> Integer
    Op = Op ⊗ eye(d)
    Op = A.vec' * Op * A.vec
    e = A.val .- A.val[1]
    res = exp.(-β*e) .* diag(Op) |> sum
    z = partitian(β, A)
    return res/z
end

function thermal_average(β::Real, Op::AbstractMatrix, A::FTLM)
    d1 = size(Op)[1]; d2 = size(A.vec)[1]
    d = d2/d1 |> Integer
    Op = Op ⊗ eye(d)
    n = size(A.vec)[1] / A.R
    res = 0.0
    for r = 1: A.R, j = 1:A.M
        e = A.val[j,r] - A.val[1,r]
        fac = (A.initv[:,r]' * A.vec[:,j,r]) * 
              (A.vec[:,j,r]' * Op * A.initv[:,r])
        res += exp(-β*e) * fac
    end
    return res*n/partitian(β, A)
end

function thermal_average(β::Real, Op::AbstractMatrix, A::OFTLM)
    d1 = size(Op)[1]; d2 = size(A.vec)[1]
    d = d2/d1 |> Integer
    Op = Op ⊗ eye(d)
    n = (size(A.vec)[1] - A.Ne) / A.R
    res = 0.0
    for r = 1: A.R, j = 1: A.M
        e = A.val[j,r] - A.eval[1]
        fac = (A.initv[:,r]' * A.vec[:,j,r]) * 
              (A.vec[:,j,r]' * Op * A.initv[:,r])
        res += exp(-β*e)*fac
    end
    res = res * n
    Op = A.evec' * Op * A.evec
    e = A.eval .- A.eval[1]
    res0 = exp.(-β*e) .* diag(Op) |> sum
    return (res+res0)/partitian(β, A)
end


"""
Correlations and Structure factor
"""
function c_average(N::Integer, m::MODEL; cf::Function, mf::Function)
    ave = 0.0
    for i = 1:N
        h = mf(m)
        ave += cf(m)
    end
    return ave/N
end

"""
The local two-time correlation functions
"""
function correlation_2time(τ::Real, β::Real,
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

"""
A useful function: f = e^(-b*e1) - e^(-b*e2) / (e2 - e1)
"""
function diffaddexp(b::Real, e1::Real, e2::Real)
    if abs(e2 - e1) < 1.e-10
        return exp(-b*e1) * b
    else
        num = exp(-b*e1) - exp(-b*e2)
        den = e2 - e1
        return num/den
    end
end


"""
Masubara frequency Green's functions: defalt type = :b
"""
function Masubara_freq_GF(n::Integer, β::Real,
    O1::T, O2::T, A::FED) where T<:AbstractMatrix
    λ = 1.0
    ωn = Masubara_freq(n,β,type=:b)
    d = d2/d1 |> Integer
    O1 = O1 ⊗ eye(d) ; O2 = O2 ⊗ eye(d)
    O1 = A.vec' * O1 * A.vec
    O2 = A.vec' * O2 * A.vec
    e = A.val .- A.val[1]
    num = 0.0
    if ωn != 0
        for i = 1: d2, j = 1: d2
            up = exp(-β*e[i]) - λ*exp(-β*e[j])
            up = up * O1[i,j] * O2[j,i]
            down = 1.0im * ωn - e[j] + e[i]
            num += up/down
        end
    else
        for i = 1: d2, j = 1: d2
            num -= O1[i,j] * O2[j,i]*diffaddexp(β,e[i],e[j])
        end
    end
    return num/partitian(β,A)
end


"""
spectral density: ρ(ω) = 2Imχ(ω) = -2ImG(ω)
"""
function spectral_density(ω::Real, β::Real, O1::T, O2::T, 
    A::FED; η::Float64 = 0.001) where T<:AbstractMatrix
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
    return 2π*num/partitian(β,A)
end


"""
structure factor(spectral representation)
"""
function structure_factor(ω::Real, β::Real,
    O1::T, O2::T, A::FED; η::Real = 0.001) where T<:AbstractMatrix
    d1 = size(O1)[1]; d2 = length(A.val)
    d = d2/d1 |> Integer
    O1 = O1 ⊗ eye(d) ; O2 = O2 ⊗ eye(d)
    O1 = A.vec' * O1 * A.vec
    O2 = A.vec' * O2 * A.vec
    e = A.val .- A.val[1]
    num = 0.0
    for i = 1: d2, j = 1: d2
        num += exp(-β*e[i])*O1[i,j]*O2[j,i]*delta(ω+e[i]-e[j],η)
    end
    return  2π*num/partitian(β,A)
end

function structure_factor(ω::Real, β::Real,
    O1::T, O2::T, A::FTLM; η::Real = 0.001) where T<:AbstractMatrix
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
        T2, Q2 = itFOLM(A.model.h, v0, nev = A.M)
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



"""
-----success, need to be modified-----
function correlation_2time(τ::Real, β::Real,
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
        T2, Q2 = itFOLM(A.model.h, v0, nev = A.M)
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


"""
-------fail!!!!---------
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
        T2, Q2 = itFOLM(A.model.h, [A.evec v0], nev = A.M-1)
        Q2 = [v0 Q2]
        T2 = Q2' * A.model.h * Q2
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

#end  #module PhysicalObservables