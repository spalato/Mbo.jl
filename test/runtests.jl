using Mbo
using Base.test


@testset "Hilbert paths" begin
    @testset "equality" begin
        p1 = HilbertPath("G", "A", "G")
        p2 = HilbertPath("G", "A", "G")
        p3 = HilbertPath("G", "B", "G")
        p4 = HilbertPath("G", "B", "G", "A", "G")
        @test p1 == p2 # despite ! p1 === p2
        @test p1 != p3
        @test p1 != p4
        @test order(p1) == 1
        @test order(p4) == 3
    end
end
#=
@testset "System" begin
    @testset "Generation" begin
        # make sure the lookups are idempotent
        tags = "GABCF"
        e_dict = Dict(zip(tags, 1:length(tags)))
        dipoles = Dict()
        for i=tags, j=tags
            if i==j dipoles[i,j] = 0
            elseif (i in "GF" && j in "GF") dipoles[i,j] = 0
            elseif (i in "ABC" && j in "ABC") dipoles[i,j] = 0
                else dipoles[i,j]=1
            end
        end
        lineshapes = Dict()
        for i=tags, j=tags
            
        end
        s = System()


    end
end
=#
#=
@testset "Correlation functions" begin
    @testset "Gridded" begin
        ts = linspace(0, 100, 101)
        dt = step(ts)
        @assert dt == 1.0
        @. cf = exp(-ts/20)+exp(-0.5*(ts/20)^2) + 0.2*cos(ts*2*pi/10)
    end
end
=#
