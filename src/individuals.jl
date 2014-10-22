using StatsBase: WeightVec, sample

immutable Female <: Organism
    auto::(Gene, Gene)
    x::(Gene, Gene)
    mito::Gene
end
Female() = Female((Gene(0 ,0), Gene(0, 0)), (Gene(0, 0), Gene(0, 0)), Gene(0, 0))

immutable Male <: Organism
    auto::(Gene, Gene)
    x::Gene
    y::Gene
end
Male() = Male((Gene(0, 0), Gene(0, 0)), Gene(0, 0), Gene(0, 0))

function Female(t, mom::Female, pop::Male, mut::Float64, rec::GeneStateRecorder)
    g = Array(Gene, 5)
    g[1] = rand() < 0.5 ? mom.auto[1] : mom.auto[2]
    g[2] = rand() < 0.5 ? mom.x[1] : mom.x[2]
    g[3] = mom.mito
    g[4] = rand() < 0.5 ? pop.auto[1] : pop.auto[2]
    g[5] = pop.x
    for i = 1:5
        g[i] = rand() < mut ? mutate!(t, g[i], rec) : g[i]
    end
    Female((g[1], g[4]), (g[2], g[5]), g[3])
end

function Male(t, mom::Female, pop::Male, mut::Float64, rec::GeneStateRecorder)
    g = Array(Gene, 4)
    g[1] = rand() < 0.5 ? mom.auto[1] : mom.auto[2]
    g[2] = rand() < 0.5 ? mom.x[1] : mom.x[2]
    g[3] = rand() < 0.5 ? pop.auto[1] : pop.auto[2]
    g[4] = pop.y
    for i = 1:4
        g[i] = rand() < mut ? mutate!(t, g[i], rec) : g[i]
    end
    Male((g[1], g[3]), g[2], g[4])
end

function learn(deme, f, m, fr::Float64, mr::Float64)
    if deme == f == m
        convert(DemeIndex, deme)
    else
        r = (deme == f ? 1.0 : 1 - fr) * (deme == m ? 1.0 : 1 - mr)
        convert(DemeIndex, rand() < r ? deme : 3 - deme)
    end
end

fitness(deme, trait, r1, r2) = convert(Float64, deme == trait ? r1 : r2)

function migrate!(t, p::Female, rec::GeneStateRecorder)
    a = (migrate!(t, p.auto[1], rec), migrate!(t, p.auto[2], rec))
    x = (migrate!(t, p.x[1], rec), migrate!(t, p.x[2], rec))
    m = migrate!(t, p.mito, rec)
    Female(a, x, m)
end

function migrate!(t, p::Male, rec::GeneStateRecorder)
    a = (migrate!(t, p.auto[1], rec), migrate!(t, p.auto[2], rec))
    x = migrate!(t, p.x, rec)
    y = migrate!(t, p.y, rec)
    Male(a, x, y)
end

function pickparent(w1::WeightVec, w2::WeightVec, r::Float64)
    convert(Int, rand() < r ? sample(w1) : sample(w2) + length(w1))
end
