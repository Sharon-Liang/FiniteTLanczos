#module FullOrthoLanczos
"""
    `randvec(N)`: Generate a N-dimesional random vector whose elements in range `1/√N (-1.0, 1.0)`
"""
function randvec(N::Integer)
    vr = map(x->x - 0.5, rand(Float64, N))
    return vr/norm(vr)
end


"""
    Iterative Classical Gram-Schmidt Algorithm.
    Input:  `u` : the vector to be orthogonalized(normalized)
            `Q` : the orthonormal basis(normalized)
    Output: `u` : the orthogonalized vector
"""
function icgs(u, Q; itmax::Integer = 3) 
    a = 0.5
    r0 = norm(u); r1 = r0
    for it = 1: itmax
        u = u - Q * (Q' * u)
        r1 = norm(u)
        if r1 > a * r0
            break
        end
        r0 = r1
    end
    if r1 <= a * r0
        @warn "Warning: Loss of Orthogonality!"
    end
    return u
end


"""
    Iterative Full Orthogonalized Lanczos Method
    Input: `A`: Symmetric Matrix
           `v0`: orthonormal basis vectors 
           `ncv`:= Number of Krylov vectors used in the computation
    Output: Tm := tridiagonal matrix
            Kb := Orthonormal basis of Krylov space
"""
function itFOLM(A, v0, lb; ncv::Integer = 50) 
    dim = size(A)[1]
    ncv = min(ncv, dim)
    Tm, Kb = zeros(ncv, ncv), zeros(dim, ncv)

    w = v0
    for k =1: ncv
        Kb[:, k] = w
        r = A * w
        Tm[k,k] = w' * r
        r = icgs(r, Kb)
        b = norm(r)
        w = r/b
        if k < ncv
            Tm[k, k+1] = b
            Tm[k+1, k] = b
        end
    end
    return Tm, Kb
end

function itFOLM(A, v0::AbstractVector; ncv::Integer = 50) 
    dim = size(A)[1]
    ncv = min(ncv, dim)
    Tm, Kb = zeros(ncv, ncv), zeros(dim, ncv)

    w = v0
    for k =1: ncv
        Kb[:, k] = w
        r = A * w
        Tm[k,k] = w' * r
        r = icgs(r, Kb)
        b = norm(r)
        w = r/b
        if k < ncv
            Tm[k, k+1] = b
            Tm[k+1, k] = b
        end
    end
    return Tm, Kb
end

function itFOLM(A ; ncv::Integer = 50) 
    v0 = randvec(size(A)[1])
    return itFOLM(A, v0, ncv = ncv) 
end



function itFOLM(A::AbstractMatrix, lb::AbstractMatrix; ncv::Integer = 50)
    """Iterative Full Orthogonalized Lanczos Method
       Input: A:= Symmetric Matrix
              nls:= number of Lanczos steps
              lb := exact low lying eigen-states
       Output: T := tridiagonal matrix
               Q := Orthonormal basis of Krylov space
    """
    dim = size(A)[1]; nls = min(nls, dim);
    ncv = min(nls, dim);
    v0 = randvec(dim) # random initiation vector
    Nv = size(lb)[2];
    v0 = icgs(v0, lb);
    Tm, Kb = zeros(ncv, ncv), zeros(dim, ncv + Nv);
    Q[:,1 : Nv] = lb
    w = v0; r = zeros(dim); k = 0;
    for k =1: ncv
        Q[:, Nv + k] = w;
        r = A * w;
        Tri[k,k] = w' * r;
        r = icgs(r, Q)
        b = norm(r)
        w = r/b;
        if k < ncv
            Tri[k, k+1] = b; Tri[k+1, k] = b
        end
    end
    #T = Q' * A * Q;
    return Tm, Kb[:, Nv+1 : ncv + Nv]
end


    
#end # module FullOrthoLanczos
