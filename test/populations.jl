srand(1234)
rates = 0.0:0.2:1.0
n = 10^4
tol = 5e-2

pop = Population([5, 5], [5, 5])
@test isa(pop[Female], SexSims.Organisms{Female})
@test isa(pop[Male], SexSims.Organisms{Male})
for s = [Female, Male]
    @test pop[s].size == (5, 5)
    @test length(pop[s].trait) == 10
    @test length(pop[s].data) == 10
    for i = 1:10
        deme = convert(SexSims.DemeIndex, i <= 5 ? 1 : 2)
        @test pop[s].trait[i] == deme
    end
end
