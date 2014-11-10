using StatsBase: WeightVec
srand(1234)
n = 10^4
rates = 0.0:0.2:1.0
tol = 5e-2

f = Female()
m = Male()
fmig = migrate!(1, f, GeneStateRecorder(5))
mmig = migrate!(1, m, GeneStateRecorder(4))

for c = [:auto, :x]
    @test all([id(g) == 0 for g = getfield(f, c)])
    @test all([mig(g) == 0 for g = getfield(f, c)])
    @test all([id(g) == 0 && mig(g) != 0 for g = getfield(fmig, c)])
end
@test id(f.mito) == 0
@test mig(f.mito) == 0
@test id(fmig.mito) == 0 && mig(fmig.mito) != 0

@test all([id(g) == 0 for g = m.auto])
@test all([mig(g) == 0 for g = m.auto])
@test all([id(g) == 0 && mig(g) != 0 for g = mmig.auto])
for c = [:x, :y]
    @test id(getfield(m, c)) == 0
    @test mig(getfield(m, c)) == 0
    @test id(getfield(mmig, c)) == 0 && mig(getfield(mmig, c)) != 0
end

for rate = rates
    fs = [Female(1, f, m, rate, GeneStateRecorder(5n)) for _ = 1:n]
    for c = [:auto, :x], i = 1:2
        frac = countnz([id(getfield(f, c)[i]) != 0 for f = fs]) / n
        @test_approx_eq_eps frac rate tol
    end
    frac = countnz([id(f.mito) != 0 for f = fs]) / n
    @test_approx_eq_eps frac rate tol

    ms = [Male(1, f, m, rate, GeneStateRecorder(4n)) for _ = 1:n]
    for i = 1:2
        frac = countnz([id(m.auto[i]) != 0 for m = ms]) / n
        @test_approx_eq_eps frac rate tol
    end
    frac = countnz([id(m.x) != 0 for m = ms]) / n
    @test_approx_eq_eps frac rate tol
    frac = countnz([id(m.y) != 0 for m = ms]) / n
    @test_approx_eq_eps frac rate tol
end

for rate = rates
    for i = 1:2, j = 1:2, k = 1:2
        data = [SexSims.learn(i, j, k, rate, 1 - rate) for _ = 1:n]
        frac = countnz(data .== i) / n
        if i == j == k
            @test frac == 1.0
        elseif i == j
            @test_approx_eq_eps frac rate tol
        elseif i == k
            @test_approx_eq_eps frac (1 - rate) tol
        else
            @test_approx_eq_eps frac (rate * (1 - rate)) tol
        end
    end
end
