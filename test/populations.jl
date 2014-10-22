srand(1234)
rates = 0.0:0.2:1.0
n = 10^4
tol = 5e-2

pop = Population((5, 5), (5, 5))
@test isa(pop[Female], SexSims.Organisms{Female})
@test isa(pop[Male], SexSims.Organisms{Male})
for s = [Female, Male]
    @test pop[s].size == (5, 5)
    @test length(pop[s].fit) == 10
    @test length(pop[s].trait) == 10
    @test length(pop[s].data) == 10
    for i = 1:10
        deme = convert(SexSims.DemeIndex, i <= 5 ? 1 : 2)
        @test pop[s].trait[i] == deme
    end
end

pids = SexSims.ParentIds(10)
@test size(pids.data) == (10, 2)
pids = SexSims.ParentIds((5, 5))
@test size(pids.data) == (10, 2)

par = Population((2, 2), (2, 2))
c = Population((2, 2), (2, 2))
pids = SexSims.ParentIds(4)
frac1, frac2 = 1.0, 1.0
data1 = SexSims.DemeIndex[0 for i = 1:2n]
data2 = SexSims.DemeIndex[0 for i = 1:2n]
for rate = rates
    pids.data[:,:] = [1 1; 1 1; 1 1; 1 1]
    for i = 1:n
        SexSims.learn!(c[Female], par, pids, (rate, 1 - rate), (1 - rate, rate))
        data1[2(i-1)+1:2i] = c[Female].trait[1:2]
        data2[2(i-1)+1:2i] = c[Female].trait[3:4]
    end
    frac1 = countnz(data1 .== 1) / 2n
    frac2 = countnz(data2 .== 1) / 2n
    @test_approx_eq_eps frac1 1.0 tol
    @test_approx_eq_eps frac2  1 - rate * (1 - rate) tol
end

p = Population((5, 5), (5, 5))
p[Female].trait[:] = [1, 1, 1, 2, 2, 1, 1, 2, 2, 2]
for rate = rates
    fitness!(p[Female], rate, 1 - rate)
    @test p[Female].fit == [1.0, 1.0, 1.0, rate, rate, 1 - rate, 1 - rate, 1.0, 1.0, 1.0]
end

pids = SexSims.ParentIds(10)
par = Population((5, 5), (5, 5))
c = Population((5, 5), (5, 5))
rec = GeneStateRecorder(1000)
par[Female].fit[:] = [1, 1, 1, 0.5, 0.5, 1.5, 1.5, 1, 1, 1]
par[Male].fit[:] = [1.5, 1.5, 1.5, 1, 1, 1, 1, 0.5, 0.5, 0.5]

dataf = [0 for _ = 1:10n]
datam = [0 for _ = 1:10n]
for rate = rates
    for i = 1:n
        SexSims.pickparents!((5, 5), pids, par, (rate, 1 - rate), (1 - rate, rate))
        SexSims.reproduce!(1, c.f, par, pids, rate, rec)
        dataf[10(i-1)+1:10i] = pids.data[:,1]
        datam[10(i-1)+1:10i] = pids.data[:,2]
        for j = 1:10
            cond = j <= 5 ? <= : >
            f = cond(pids[j, 1], 5) ? == : !=
            @test f(mig(c[Female].data[j].auto[1]), 0x0)
            @test f(mig(c[Female].data[j].x[1]), 0x0)
            @test f(mig(c[Female].data[j].mito), 0x0)
            f = cond(pids[j,2], 5) ? == : !=
            @test f(mig(c[Female].data[j].auto[2]), 0x0)
            @test f(mig(c[Female].data[j].x[2]), 0x0)
        end
    end
    fracf = countnz(dataf .<= 5) / 10n
    @test_approx_eq_eps fracf rate tol
    fracm = countnz(datam .<= 5) / 10n
    @test_approx_eq_eps fracm 1 - rate tol

    fracf = countnz(dataf .<= 3) / 10n
    predf = 3rate/4
    @test_approx_eq_eps fracf predf tol
    fracm = countnz(datam .<= 3) / 10n
    predm = 9(1-rate)/13
    @test_approx_eq_eps fracm predm tol

    fracf = countnz(dataf .>= 8) / 10n
    predf = (1-rate)/2
    @test_approx_eq_eps fracf predf tol
    fracm = countnz(datam .>= 8) / 10n
    predm = 3rate/7
    @test_approx_eq_eps fracm predm tol
end
