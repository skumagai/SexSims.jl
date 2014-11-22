using SexSims: learn, nmigrants!, nbefore!

srand(1)
n = 10^4
rates = 0.0:0.2:1.0
tol = 5e-2

for rate in rates
    for i = 1:2, j = 1:2, k = 1:2
        data = [learn(i, j, k, rate, 1 - rate) for _ = 1:n]
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

for rate in rates
    data = Array(Float64, n, 2)
    mig = Array(Int, 2)
    pop = 10000
    for i = 1:n
        nmigrants!(mig, [pop, pop], [rate, 1 - rate])
        data[i,:] = mig ./ pop
    end
    @test_approx_eq_eps mean(data[:,1]) rate tol
    @test_approx_eq_eps mean(data[:,2]) 1 - rate tol
end

nb = Array(Int, 2)
pop = 10000
nbefore!(nb, [pop, 2pop], [100, 200])
@test nb == [pop - 100, 2pop + 100]
