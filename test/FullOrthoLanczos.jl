using Test, FiniteTLanczos
using Random; Random.seed!()

@testset "icgs" begin
    N = 1024
    ortho = Vector{Bool}(undef, N)
    Q = rand(N); Q = Q/norm(Q)
    for i = 1:N
        v = rand(N); v = v/norm(v)
        u = icgs(v, Q)
        r = Q' * u
        ortho[i] = map(x->isapprox(x, 0, atol = 1.e-12), r) |> all 
        Q = hcat(Q, u/norm(u))
    end
    @test all(ortho)
end

