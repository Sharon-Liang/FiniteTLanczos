#module utilities
import Base: kron

"""
Pauli matrices
"""
function pauli(symbol::Symbol)
    if symbol==:x return [0. 1.; 1. 0.]
    elseif symbol==:y return [0. -1im; 1im 0.]
    elseif symbol==:iy return [0. 1.; -1. 0.]
    elseif symbol==:z return [1. 0.; 0. -1.]
    elseif symbol==:+ return [0. 1.; 0. 0.]
    elseif symbol==:- return [0. 0.; 1. 0.]
    else
        error("The input should be :x,:y,:z,:+,:-, :iy.")
    end
end

"""
identity matrix
"""
function eye(dim::Integer)
    Matrix{Float64}(I,dim,dim)
end

function eye(T::DataType, dim::Integer)
    Matrix{T}(I,dim,dim)
end

function ⊗(A::AbstractArray, B::AbstractArray)
    kron(A,B)
end

"""
Gaussian approximate delta function
"""
function delta(x::Real, η::Real)
    num = η
    den = (x^2 + η^2) * π
    return num/den
end

"""
Create a N-dimensional normalized random vector
"""
function random_init(dim::Integer)
    Random.seed!()
    vr = rand(dim) .- 0.5
    vr = vr/norm(vr)
    return vr
end

function random_init(T::DataType, dim::Integer)
    Random.seed!()
    vr = rand(T, dim) .- 0.5
    vr = vr/norm(vr)
    return vr
end

"""
Average deviation form full sampling.
    dim: dimension of the target vector
    R: number of random vectors
"""
function fidelity(dim::Integer, R::Integer;)
    v = zeros(dim)
    for r = 1:R
        v0 = random_init(dim)
        v += v0 .* v0
    end
    return abs.(v ./ R .- ones(dim) ./ dim ) |> sum
end

function fidelity(T::DataType, dim::Integer, R::Integer)
    v = zeros(T, dim)
    for r = 1:R
        v0 = random_init(T, dim)
        v += v0 .* v0
    end
    return abs.(v ./ R .- ones(T,dim) ./ dim ) |> sum
end


"""
Iterative Classical Gram-Schmidt Algorithm.
    Input: u := the vector to be orthogonalized
        Q := the orthonormal basis
    Output: u := the orthogonalized vector
"""
function icgs(u::AbstractArray, Q::AbstractArray; itmax::Integer = 3)
    a = 0.5
    r0 = norm(u); r1 = r0
    for it = 1: itmax
        u = u - Q * (Q' * u)
        r1 = norm(u)
        if r1 > a * r0
            break
        end
        it += 1; r0 = r1
    end
    if r1 <= a * r0
        @warn "Warning: Loss of Orthogonality!"
    end
    return u
end


"""
Iterative Full Orthogonalized Lanczos Method
    Input: A:= Symmetric Matrix
        nev:= number of Lanczos steps
    Output: T := tridiagonal matrix
            Q := Orthonormal basis of Krylov space
"""
function itFOLM(A::AbstractMatrix; nev::Integer = 50, return_basis::Bool=true)
    dim = size(A)[1]; nev = min(nev, dim);
    ncv = min(nev, dim);
    v0 = random_init(dim) # random initiation vector
    T, Q = zeros(ncv, ncv), zeros(dim, ncv);
    w = v0; r = zeros(dim); k = 0;
    for k =1: ncv
        Q[:, k] = w;
        r = A * w;
        T[k,k] = w' * r;
        r = icgs(r, Q)
        b = norm(r)
        w = r/b;
        if k < ncv
            T[k, k+1] = b; T[k+1, k] = b
        end
    end
    #T = Q' * A * Q;
    return_basis ? (return T,Q) : (return T)
end


"""
Iterative Full Orthogonalized Lanczos Method
    Input: A:= Symmetric Matrix
           nev:= number of Lanczos steps
           v0: init vector
    Output: T := tridiagonal matrix
            Q := Orthonormal basis of Krylov space
"""
function itFOLM(A::AbstractMatrix, v0::AbstractVector; nev::Integer = 50, return_basis::Bool=true)
    dim = size(A)[1]; nev = min(nev, dim);
    ncv = min(nev, dim);
    T, Q = zeros(ncv, ncv), zeros(dim, ncv);
    w = v0; r = zeros(dim); k = 0;
    for k =1: ncv
        Q[:, k] = w;
        r = A * w;
        T[k,k] = w' * r;
        r = icgs(r, Q)
        b = norm(r)
        w = r/b;
        if k < ncv
            T[k, k+1] = b; T[k+1, k] = b
        end
    end
    #T = Q' * A * Q;
    return_basis ? (return T,Q) : (return T)
end


"""
Iterative Full Orthogonalized Lanczos Method
    Input: A:= Symmetric Matrix
        nev:= number of Lanczos steps
        lb := exact low lying eigen-states
    Output: T := tridiagonal matrix
            Q := Orthonormal basis of Krylov space
"""
function itFOLM(A::AbstractMatrix, lb::AbstractMatrix; nev::Integer = 50, return_basis::Bool=true)
    dim = size(A)[1]; nev = min(nev, dim);
    ncv = min(nev, dim);
    v0 = random_init(dim) # random initiation vector
    Nv = size(lb)[2];
    v0 = icgs(v0, lb);
    T, Q = zeros(ncv, ncv), zeros(dim, ncv + Nv);
    Q[:,1 : Nv] = lb
    w = v0; r = zeros(dim); k = 0;
    for k =1: ncv
        Q[:, Nv + k] = w;
        r = A * w;
        T[k,k] = w' * r;
        r = icgs(r, Q)
        b = norm(r)
        w = r/b;
        if k < ncv
            T[k, k+1] = b; T[k+1, k] = b
        end
    end
    #T = Q' * A * Q;
    return_basis ? (return T, Q[:, Nv+1 : ncv + Nv]) : (return T)
end


#end module utilities
