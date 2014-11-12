srand(1234)
rates = 0.0:0.2:1.0
n = 10^4
tol = 5e-2

pop = Population([5, 5], [5, 5])
@test isa(pop[Female], SexSims.Organisms{Female})
@test isa(pop[Male], SexSims.Organisms{Male})
for s = [Female, Male]
    @test pop[s].size == [5, 5]
    @test length(pop[s].trait) == 10
    @test length(pop[s].data) == 10
    for i = 1:10
        deme = convert(SexSims.DemeIndex, i <= 5 ? 1 : 2)
        @test pop[s].trait[i] == deme
    end
end

for i = 0:20, j = 0:i, _ = 1:1000
    d = 20 - i + j

    nelem = [i, 40 - i]
    sels = SexSims.selectmigrants(nelem, [j, d])
    @test length(sels[1]) == length(sels[2]) - 20 + i
    @test all(sels[1] .<= i)
    @test all(sels[2] .> i)
    elems = [1:40]
    tar, src = SexSims.listmovers(nelem, [j, d])
    elems[tar] = elems[src]
    @test length(unique(elems)) == 40
    @test countnz((elems .> i)[1:20]) == d
    @test countnz((elems .<= i)[21:end]) == j
    rec = GeneStateRecorder(100)
    didx = vcat(fill(convert(SexSims.DemeIndex, 1), i), fill(convert(SexSims.DemeIndex, 2), 40 - i))
    orgs = SexSims.Organisms{Female}(nelem, didx, [Female() for _ = 1:sum(nelem)])
    migrate!(1, orgs, nelem, [j, d], rec)
    @test countnz([SexSims.mig(o.auto[1]) != 0x0 for o in sub(orgs.data, 1:20)]) == d
    @test countnz([SexSims.mig(o.auto[1]) != 0x0 for o in sub(orgs.data, 21:40)]) == j

    nelem = [40 - i, i]
    sels = SexSims.selectmigrants(nelem, [d, j])
    @test length(sels[1]) == length(sels[2]) + 20 - i
    @test all(sels[1] .<= 40 - i)
    @test all(sels[2] .> 40 - i)
    elems = [1:40]
    tar, src = SexSims.listmovers(nelem, [d, j])
    elems[tar] = elems[src]
    @test length(unique(elems)) == 40
    @test countnz((elems .> 40 - i)[1:20]) == j
    @test countnz((elems .<= 40 - i)[21:end]) == d
    rec = GeneStateRecorder(100)
    didx = vcat(fill(convert(SexSims.DemeIndex, 1), 40 - i), fill(convert(SexSims.DemeIndex, 2), i))
    orgs = SexSims.Organisms{Female}(nelem, didx, [Female() for _ = 1:sum(nelem)])
    migrate!(1, orgs, nelem, [d, j], rec)
    @test countnz([SexSims.mig(o.auto[1]) != 0x0 for o in sub(orgs.data, 1:20)]) == j
    @test countnz([SexSims.mig(o.auto[1]) != 0x0 for o in sub(orgs.data, 21:40)]) == d
end
