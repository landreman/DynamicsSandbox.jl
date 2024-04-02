using Test
using DynamicsSandbox

@testset "test derivative of standard map" begin
    x = 0.3
    y = -0.2
    k = 0.9
    D = DynamicsSandbox.standard_map_gradient(x, y, k)

    ϵ = 1e-8
    vp = DynamicsSandbox.standard_map(x + ϵ, y, k)
    vm = DynamicsSandbox.standard_map(x - ϵ, y, k)
    @test (vp - vm) / (2ϵ) ≈ D[:, 1]
    vp = DynamicsSandbox.standard_map(x, y + ϵ, k)
    vm = DynamicsSandbox.standard_map(x, y - ϵ, k)
    @test (vp - vm) / (2ϵ) ≈ D[:, 2]
end

@testset "test inverse of standard map" begin
    k = 0.9
    x = 0.3
    y = 0.2
    xp, yp = DynamicsSandbox.standard_map(x, y, k)
    xpp, ypp = DynamicsSandbox.standard_map_inverse(xp, yp, k)
    @test x ≈ xpp
    @test y ≈ ypp
end
