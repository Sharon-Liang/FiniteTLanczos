#module methods
abstract type ED end

struct FED <: ED
    m::model
    val::AbstractVector
    vec::AbstractMatrix
end

function FED(A::model)
    H = Array(A.ham)
    e, v = eigen(H)
    return FED(A, e, v)
end

struct FTLM <: ED
    m::model
    R::Integer
    M::Integer
    initv::AbstractArray
    val::AbstractMatrix
    vec::AbstractArray
end

struct OFTLM <:ED
    m::model
    Ne::Integer # number of exact low lying eigenstates
    eval::AbstractVector
    evec::AbstractMatrix
    R::Integer
    M::Integer
    initv::AbstractArray
    val::AbstractMatrix
    vec::AbstractArray
end

function FTLM(A::model; R::Integer = 50, M::Integer = 90)
    """Finite Temperature Lanczos Method
       Input: A := Hamiltonian Matrix (Hermitian)
              M := The number of Lanczos step
              R := The number of random sampling
        Output: initv[dim,R], val[M,R], vec[dim,M,R] |> FTLM
    """
    dim = size(A.ham)[1]
    initv = zeros(dim, R)
    ncv = min(M, dim)
    val = zeros(ncv,R); vec = zeros(dim, ncv, R)
    for r = 1:R
        T, Q = itFOLM(A.ham, nev = ncv)
        initv[:,r] = Q[:,1]
        e, v = eigen(T)
        val[:,r] = e; vec[:,:,r] = Q * v
    end
    return FTLM(A, R, ncv, initv, val, vec)
end

function OFTLM(A::model; R = 50, M = 90, Ne =10)
    """Orthogonalized Finite Temperature Lanczos Method
       Input: A := Hamiltonian Matrix
              M := The number of Lanczos step
              R := The number of random sampling
              Ne := Number of Exact eigenstates
        Output: V := [E(rj),  <v psi>*<psi v>, <v psi>*<psi O v>]
                dim/R
    """
    dim = size(A.ham)[1]; 
    initv = zeros(dim, R)
    ncv = min(M, dim)
    val = zeros(ncv,R); vec = zeros(dim, ncv, R)

    Ne = min(Ne, dim-1)
    Ee, Ve = eigs(A.ham, nev = Ne, which =:SR)
    
    for r = 1: R
        T, Q = itFOLM(A.ham, Ve, nev = ncv)
        initv[:,r] = Q[:,1]
        e, v = eigen(T)
        val[:,r] = e; vec[:,:,r] = Q * v
    end 
    return OFTLM(A, Ne, Ee, Ve, R, ncv, initv, val, vec)
end

#end  # module methods
