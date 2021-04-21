#module utilities
import Base: kron

function pauli(char::Char)
    if char=='x' return [0. 1.; 1. 0.] |> sparse
    elseif char=='y' return [0. -1im; 1im 0.] |> sparse
    elseif char=='z' return [1. 0.; 0. -1.] |> sparse
    elseif char=='+' return [0. 1.; 0. 0.] |> sparse
    elseif char=='-' return [0. 0.; 1. 0.] |> sparse
    else
        error("The input should be 'x','y','z','+','_'.")
    end
end

function eye(N::Integer)
    Matrix(1.0I, N, N) |> sparse
end

function ⊗(A::AbstractMatrix, B::AbstractMatrix)
    kron(A,B)
end

function random_init(N::Integer)
    vr = rand(Float64, N) .- 0.5
    vr = vr/norm(vr)
    return vr
end

function icgs(u::AbstractArray, Q::AbstractArray)
    """Iterative Classical Gram-Schmidt Algorithm.
       Input: u := the vector to be orthogonalized
              Q := the orthonormal basis
       Output: u := the orthogonalized vector
    """
    a = 0.5; itmax = 3;
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
        println("Warning: Loss of Orthogonality!")
    end
    return u
end

function itFOLM(A::AbstractMatrix; nev::Integer = 50, lb = nothing, return_basis = true)
    """Iterative Full Orthogonalized Lanczos Method
       Input: A:= Symmetric Matrix
              nev:= number of Lanczos steps
              lb := exact low lying eigen-states
       Output: T := tridiagonal matrix
               Q := Orthonormal basis of Krylov space
    """
    dim = size(A)[1]; nev = min(nev, dim);
    ncv = min(nev, dim);
    v0 = random_init(dim) # random initiation vector
    if lb == nothing
        Nv = 0
        T, Q = zeros(ncv, ncv), zeros(dim, ncv + Nv);
    else
        Nv = size(lb)[2];
        v0 = icgs(v0, lb);
        T, Q = zeros(ncv, ncv), zeros(dim, ncv + Nv);
        Q[:,1 : Nv] = lb
    end


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
    if return_basis
        return T, Q[:, Nv+1 : ncv + Nv]
    else
        return T
    end
end
#end  # module utilities
