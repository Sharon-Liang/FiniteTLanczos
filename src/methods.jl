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

#struct FTLM{T<:Integer}
#    R::T
#    M::T
#    initv::AbstractMatrix
#    val::AbstractMatrix
#    vec::AbstractArray{T1,3} where T1<:Number
#end



#function FTLM(A::AbstractMatrix; R::Integer = 50, M::Integer = 90)
#    """Partition Function by the Finite Temperature Lanczos Method
#       Input: A := Hamiltonian Matrix (Hermitian)
#              M := The number of Lanczos step
#              R := The number of random sampling
#        Output: initv[dim,M,R], val[M,R], vec[M,M,R] |> FTLM
#    """
#    dim = size(A)[1]
#    initv = zeros(dim,R)
#    val = zeros(M,R); vec = zeros(ComplexF64, dim, M, R)
#    for r = 1:R
#        T, Q = itFOLM(A, nev = M)
#        e, v = eigen(T)
#        initv[:,r] = Q[:,1]
#        val[:,r] = e; vec[:,:,r] = v
#    end
#    return FTLM(R, M, initv, val, vec)
#end

#function FTLM(A::AbstractMatrix; R::Integer = 50, M::Integer = 90, Op = nothing)
#    """Partition Function by the Finite Temperature Lanczos Method
#       Input: A := Hamiltonian Matrix
#              M := The number of Lanczos step
#              R := The number of random sampling
#              temp := Temperature
#              Op := A general operator
#        Output: V := [E(rj),  <v psi>*<psi v>, <v psi>*<psi O v>]
#                dim/R
#    """
#    dim = size(A)[1]; fac = dim/R
#    if Op == nothing
#        n = 2
#    else
#        n = 3
#    end
#
#    V = zeros(R, M, n)
#    for r = 1:R
#        T, Q = itFOLM(A, nev = M)
#        vals, vecs = eigen(T)
#        emin = minimum(vals)
#        for j = 1:M
#            V[r,j,1] = vals[j] - emin
#            V[r,j,2] = vecs[1,j] * vecs[1,j]' * fac
#            if Op != nothing
#                V[r,j,3] = vecs[1,j] * (vecs[:,j]' * Q' * Op * Q[:, 1])
#            end
#        end
#    end
#    return V
#end
#end  # module methods
