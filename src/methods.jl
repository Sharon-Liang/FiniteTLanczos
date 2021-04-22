#module methods
struct FED
    L::Integer
    val::AbstractVector
    vec::AbstractMatrix
end

function FED(A::model)
    H = Array(A.ham)
    e, v = eigen(H)
    return FED(A.L, e, v)
end

struct FTLM{T<:Integer}
    R::T
    M::T
    initv::AbstractArray
    val::AbstractMatrix
    vec::AbstractArray
end

struct OFTLM{T<:Integer}
    Ne::T # number of exact low lying eigenstates
    eval::AbstractVector
    evec::AbstractMatrix
    R::T
    M::T
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
    val = zeros(M,R); vec = zeros(dim, M, R)
    for r = 1:R
        T, Q = itFOLM(A.ham, nev = M)
        initv[:,r] = Q[:,1]
        e, v = eigen(T)
        val[:,r] = e; vec[:,:,r] = Q * v
    end
    return FTLM(R, M, initv, val, vec)
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
    Ee, Ve = eigs(A.ham, nev = Ne, which =:SR)
    dim = size(A.ham)[1]; 
    initv = zeros(dim, R)
    val = zeros(M,R); vec = zeros(dim, M, R)   
    for r = 1: R
        T, Q = itFOLM(A.ham, Ve, nev = M)
        initv[:,r] = Q[:,1]
        e, v = eigen(T)
        val[:,r] = e; vec[:,:,r] = Q * v
    end 
    return OFTLM(Ne, Ee, Ve, R, M, initv, val, vec)
end

#end  # module methods
