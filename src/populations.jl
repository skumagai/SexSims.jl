immutable Parents
    moms::(Vector{(DemeIndex, Int)}, Vector{DemeIndex, Int})
    pops::(Vector{(DemeIndex, Int)}, Vector{DemeIndex, Int})
end
function Parents(n::(Int, Int))
    f = (Array((DemeIndex, Int), n[1]), Array((DemeIndex, Int), n[2]))
    m = (Array((DemeIndex, Int), n[1]), Array((DemeIndex, Int), n[2]))
    Parents(f, m)
end
immutable ParentsIds
    tofem::Parents
    tomale::Parents
end
function ParentsIds(nf::(Int, Int), nm::(Int, Int))
    ParentsIds(Parents(nf), Parents(nm))
end

function learn(deme::DemeIndex, f::DemeIndex, m::DemeIndex, fr::Float64, mr::Float64)
    if deme == f == m
        deme
    else
        r = (deme == f ? 1.0 : 1 - fr) * (deme == m ? 1.0 : 1 - mr)
        rand() < r ? deme : 3 - deme
    end
end

function learn!{T<:Organisms}(c::T, par::Population,
    pids::ParentsIds,
    fr::(Float64, Float64),
    mr::(Float64, Float64))
    for deme = 1:2
        ct = c.trait[deme]
        momt = par.f.trait
        popt = par.m.trait
        momids = pids.tofem[deme]
        popids = pids.tomale[deme]
        for i = 1:length(ct)
            mdeme, mid = momoids[i]
            pdeme, pid = popids[i]
            ct[i] = learn(deme, momt[mdeme][mid], popt[pdeme][pid], fr[deme], mr[deme])
        end
    end
    nothing
end

function learn!(c::Population, p::Population,
        pids::ParentsIds,
        fr::(Float64, Float64),
        mr::(Float64, Float64))
    learn!(c.f, p, pids.tofem, fr, fm)
    learn!(c.m, p, pids.tomale, fr, fm)
    nothing
end

function fitness(deme::DemeIndex, trait::DemeIndex, r1::Float64, r2::Float64)
    deme == trait ? r1 : r2
end

function fitness!{T<:Organisms}(p::T, r1::Float64, r2::Float64)
    for deme = 1:2
        f = p.fit[deme]
        t = p.trait[deme]
        rd1, rd2 = deme == 1 ? (1.0, r1) : (r2, 1.0)
        for i = 1:p.size[1]
            f[i] = fitness(deme, t[i], rd1, rd2)
        end
    end
    nothing
end

function reproduction!{T<:Organisms}(c::T, mom::Females, pop::Males,
    f::Vector{(DemeIndex, Int), (DemeIndex, Int)},
    m::Vector{(DemeIndex, Int), (DemeIndex, Int)})
    for deme = 1:2
end

abstract Organism

immutable Female <: Organism
    auto::(Gene, Gene)
    x::(Gene, Gene)
    mito::Gene
end

function Female(rec::GeneStateRecorder)
    Female((getnewgene(rec), getnewgene(rec)),
        (getnewgene(rec), getnewgene(rec)),
        getnewgene(rec))
end

function Female(t::Generation, mom::Female, pop::Male, mut::Float64, rec::GeneStateRecorder)
    g = Arrray(Gene, 5)
    g[1] = rand() < 0.5 ? mom.auto[1] : mom.auto[2]
    g[2] = rand() < 0.5 ? mom.x[1] : mom.x[2]
    g[3] = mom.mito
    g[4] = rand() < 0.5 ? pop.auto[1] : pop.auto[2]
    g[5] = pop.x
    for i = 1:5
        g[i] = rand() < mut ? mutate!(t, g[i], rec) : g[i]
    end
    Female((g[1], g[4]), (g[2], g[4]), g[5])
end

immutable Male <: Organism
    auto::(Gene, Gene)
    x::Gene
    y::Gene
end

function Male(rec::GeneStateRecorder)
    Male((getnewgene(rec), getnewgene(rec)),
        getnewgene(rec),
        getnewgene(rec))
end

function Male(t::Generation, mom::Female, pop::Male, mut::Float64, rec::GeneStateRecorder)
    g = Array(Gene, 4)
    g[1] = rand() < 0.5 ? mom.auto[1] : mom.auto[2]
    g[2] = rand() < 0.5 ? mom.x[1] : mom.x[2]
    g[3] = rand() < 0.5 ? pop.auto[1] : pop.auto[2]
    g[4] = pop.y
    for i = 1:4
        g[i] = rand() < mut ? mutate!(t, g[i], rec) : g[i]
    end
    Male((g[1], g[3]), (g[2], g[4]))
end

function migrate!(t::Generation, p::Female, rec::GeneStateRecorder)
    a = (migrate!(t, p.auto[1], rec), migrate!(t, p.auto[2], rec))
    x = (migrate!(t, p.x[1], rec), migrate!(t, p.x[2], rec))
    m = migrate!(t, p.mito, rec)
    Female(a, x, m)
end

function migrate!(t::Generation, p::Male, rec::GeneStateRecorder)
    a = (migrate!(t, p.auto[1], rec), migrate!(t, p.auto[2], rec))
    x = migrate!(t, p.x, rec)
    y = migrate!(t, p.y, rec)
    Male(a, x, y)
end

abstract Organisms
type Females <: Organisms
    size::(Int, Int)
    fit::(Vector{Float64}, Vector{Float64})
    trait::(Vector{DemeIndex}, Vector{DemeIndex})
    data::Vector{Female}
end

type Males <: Organisms
    size::(Int, Int)
    fit::(Vector{Float64}, Vector{Float64})
    trait(Vector{DemeIndex}, Vector{DemeIndex})
    data::Vector{Male}
end

immutable Population
    f::Females
    m::Males
end

function Population(nf::(Int, Int), nm::(Int, Int))
    ffit = tuple([Array(Float64, n) for n = nf]...)
    mfit = tuple([Array(Float64, n) for n = nm]...)

    ft = tuple([fills(convert(DemeIndex, i), j) for (i, j) = enumerate(nf)]...)
    mt = tuple([fills(convert(DemeIndex, i), j) for (i, j) = enumerate(nm)]...)

    f = Females(nf, ffit, ft)
    m = Males(nm, mfit, mt)

    Population(f, m)
end
