using Distributions: Binomial

function learn(deme, f, m, fr::Float64, mr::Float64)
    if deme == f == m
        deme
    else
        r = (deme == f ? 1.0 : 1 - fr) * (deme == m ? 1.0 : 1 - mr)
        rand() < r ? deme : 0x3 - deme
    end
end

function nmigrants!(mig, n, v)
    for i = 1:length(n)
        mig[i] = rand(Binomial(n[i], v[i]))
    end
end

function nbefore!(nb, n, m)
    for i = 1:length(n)
        nb[i] = n[i] + m[i] - m[3 - i]
    end
end
