#module methods
abstract type AbstractLanczos end

"""
FED: full exact diagonalization
     model: physical models
     val: eigen values
     vec: eigen vectors
"""
struct FED{T<:Number} <: LANCZOS
    model::MODEL
    val::Vector{T}
    vec::Matrix{T}
end

function FED(A::MODEL)
    H = Array(A.h)
    e, v = eigen(H)
    return FED(A, e, v)
end


"""
FTLM: finite-temperature lanczos method
     model: physical models
     R: number of random samplings
     M: number of Lanczos steps
     initv: init vector
     val: eigen values
     vec: eigen vectors
"""
struct FTLM{T<:Number,Ti<:Integer} <: LANCZOS
    model::MODEL
    R::Ti
    M::Ti
    initv::Array{T}
    val::Matrix{T}
    vec::Array{T}
end

"""
FTLM: finite-temperature lanczos method
    Input: A := Hamiltonian Matrix (Hermitian)
           M := number of Lanczos steps
           R := number of random samplings
    Output: initv[dim,R], val[M,R], vec[dim,M,R] |> FTLM
"""
function FTLM(A::MODEL; R::Ti=50, M::Ti=90) where Ti<:Integer
    dim = size(A.h)[1]
    initv = zeros(dim, R)
    ncv = min(M, dim)
    val = zeros(ncv,R); vec = zeros(dim, ncv, R)
    for r = 1:R
        T, Q = itFOLM(A.h, nev = ncv)
        initv[:,r] = Q[:,1]
        e, v = eigen(T)
        val[:,r] = e; vec[:,:,r] = Q * v
    end
    return FTLM(A, R, ncv, initv, val, vec)
end


"""
OFTLM: orthogonalized finite-temperature lanczos method
     model: physical models
     Ne:number of exact low lying eigenstates
     eval:exact eigen values
     evec:exact eigen vectors
     R: number of random samplings
     M: number of Lanczos steps
     initv: init vector
     val: eigen values
     vec: eigen vectors
"""
struct OFTLM{T<:Number,Ti<:Integer} <: LANCZOS
    model::MODEL
    Ne::Ti
    eval::Vector{T}
    evec::Matrix{T}
    R::Ti
    M::Ti
    initv::Array{T}
    val::Matrix{T}
    vec::Array{T}
end

"""
Orthogonalized Finite Temperature Lanczos Method
    Input: A := Hamiltonian Matrix
        M := The number of Lanczos step
        R := The number of random sampling
        Ne := Number of Exact eigenstates
    Output: V := [E(rj),  <v psi>*<psi v>, <v psi>*<psi O v>]
            dim/R
"""
function OFTLM(A::MODEL; R::Ti=50, M::Ti=90, Ne::Ti=10) where Ti<:Integer
    dim = size(A.h)[1]; 
    initv = zeros(dim, R)
    ncv = min(M, dim)
    val = zeros(ncv,R); vec = zeros(dim, ncv, R)

    Ne = min(Ne, dim-1)
    Ee, Ve = eigs(A.h, nev = Ne, which =:SR)
    
    for r = 1: R
        T, Q = itFOLM(A.h, Ve, nev = ncv)
        initv[:,r] = Q[:,1]
        e, v = eigen(T)
        val[:,r] = e; vec[:,:,r] = Q * v
    end 
    return OFTLM(A, Ne, Ee, Ve, R, ncv, initv, val, vec)
end

#end #module methods
