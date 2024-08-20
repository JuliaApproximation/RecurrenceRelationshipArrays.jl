using RecurrenceRelationshipArrays, RecurrenceRelationships, FillArrays, LazyArrays, InfiniteArrays, BandedMatrices, Test

rec_T = Vcat(1, Fill(2,∞)), Zeros{Int}(∞), Ones{Int}(∞)
rec_U = Fill(2,∞), Zeros{Int}(∞), Ones{Int}(∞)

U = (N,z) -> forwardrecurrence(N, rec_U..., z)
T = (N,z) -> forwardrecurrence(N, rec_T..., z)

@testset "Clenshaw" begin
    X = SymTridiagonal(Fill(0.0,∞), Fill(1/2,∞))
    C = Clenshaw([1,2,3], rec_U..., X)
    @test C[1:4,1:4] isa BandedMatrix
end

@testset "RecurrenceArray" begin
    @testset "RecurrenceVector" begin()
        for z in (0.1, 1.0)
            r = RecurrenceArray(z, rec_T, [0.0, 1.0])
            @test r[1:1000] ≈ [0; U(999, z)]
            @test r[10_000] ≈ U(9_999, z)[end]

            r = RecurrenceArray(z, rec_U, [1.0, z])
            @test r[1:1000] ≈ T(1000, z)
            @test r[10_000] ≈ T(10_000, z)[end]
        end

        for z in (1.000000001, 1.000001, -1.000001, 10.0, -10.0)
            ξ = inv(z + sign(z)sqrt(z^2-1))
            r = RecurrenceArray(z, rec_U, [ξ,ξ^2])
            @test r[1:1000] ≈ ξ.^(1:1000)
            @test r[10_000] ≈ ξ.^(10_000) atol=3E-10
        end

        for z in (0.2567881003580743 - 0.33437737333561895im)
            ξ = π * (z - √(z-1) * √(z+1))
            r = RecurrenceArray(z, rec_U, [ξ,2z*ξ-π])
            @test sizeof(r.data) < 1000_000
        end
    end

    @testset "RecurrenceMatrix" begin
        z = [2.,3.,100.]
        ξ = @. inv(z + sign(z)sqrt(z^2-1))
        r = RecurrenceArray(z, rec_U, [ξ'; ξ'.^2])
        @test r[1:100,:] ≈ [RecurrenceArray(z[j], rec_U,  [ξ[j],ξ[j]^2])[k] for k=1:100, j=axes(z,1)]
        @test r[1:100,:] ≈ r[1:100,1:3] ≈ r[collect(1:100),1:3] ≈ r[1:100,collect(1:3)] ≈ r[collect(1:100),:]
        @test r[100,:] ≈ r[100,1:3]
    end
end
