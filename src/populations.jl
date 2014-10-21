immutable Organisms{T<:Organism}
    size::(Int, Int)
    fit::Vector{Float64}
    trait::Vector{DemeIndex}
    data::Vector{T}
end

immutable Population
    f::Organisms{Female}
    m::Organisms{Male}
end

function Population(nf::(Int, Int), nm::(Int, Int))
    ffit = Array(Float64, sum(nf))
    mfit = Array(Float64, sum(nm))

    ft = vcat([fills(convert(DemeIndex, i), j) for (i, j) = enumerate(nf)]...)
    mt = vcat([fills(convert(DemeIndex, i), j) for (i, j) = enumerate(nm)]...)

    f = Organisms{Female}(nf, ffit, ft, [Female() for _ = sum(nf)])
    m = Organisms{Male}(nm, mfit, mt, [Male() for _ = sum(nm)])

    Population(f, m)
end

Base.getindex(p::Population, ::Type{Female}) = p.f
Base.getindex(p::Population, ::Type{Male}) = p.m

immutable ParentIds
    data::Matrix{Int}
end
ParentIds(n::(Int, Int)) = ParentIds(Array(Int, sum(n), 2))
ParentIds(n::Int) = ParentIds(Array(Int, n, 2))

Base.getindex(p::ParentIds, i) = p.data[i, :]
Base.getindex(p::ParentIds, i, j) = p.data[i, j]

function learn!{T<:Organism}(c::Organisms{T}, par::Population,
    pids::ParentIds,
    fr::(Float64, Float64),
    mr::(Float64, Float64))
    for i = 1:length(c.data)
        deme = i <= c.size[1] ? 1 : 2
        fid, mid = pids[i]
        c.trait[i] = learn(deme, par.f.trait[fid], par.m.trait[mid], fr[deme], mr[deme])
    end
    nothing
end

function learn!(c::Population, p::Population,
    pids::(ParentIds, ParentIds),
    ftof::(Float64, Float64),
    mtof::(Float64, Float64),
    ftom::(Float64, Float64),
    mtom::(Float64, Float64))
    learn!(c.f, p, pids[1], ftof, mtof)
    learn!(c.m, p, pids[2], ftom, mtom)
    nothing
end

function fitness!{T<:Organism}(p::Organisms{T}, r1::Float64, r2::Float64)
    for i = 1:length(p.data)
        deme = i <= p.size[1] ? 1 : 2
        rd1, rd2 = deme == 1 ? (1.0, r1) : (r2, 1.0)
        p.fit[i] = fitness(deme, p.trait[i], rd1, rd2)
    end
    nothing
end

function fitness!(p::Population, fr::(Float64, Float64), mr::(Float64, Float64))
    fitness!(p.f, fr...)
    fitness!(p.m, mr...)
end

function reproduce!{T<:Organism}(t::Generation, c::Organisms{T}, p::Population,
    pid::ParentIds,
    mut::Float64,
    rec::GeneStateRecorder)

    sep = c.size[1]
    for i = 1:length(c.data)
        cdeme = i < sep ? 1 : 2
        fid, mid = pid[i]
        fdeme = fid <= sep ? 1 : 2
        mdeme = mid <= sep ? 1 : 2
        f = p.f.data[fid]
        m = p.m.data[mid]
        mom = cdeme == fdeme ? f : migrate!(t, f, rec)
        pop = cdeme == mdeme ? m : migrate!(t, m, rec)
        c[i] = T(t, mom, pop, mut, rec)
    end
end

function reproduce!(t::Generation, c::Population, p::Population,
    pid::(ParentIds, ParentIds),
    mut::Float64,
    rec::GeneStateRecorder)
    reproduce!(t, c.f, p, pid[1], mut, rec)
    reproduce!(t, c.m, p, pid[2], mut, rec)
    nothing
end

function pickparents!{T<:Organism}(::Type{T}, ids::ParentIds, pop::Population, fb::(Float64, Float64), mb::(Float64, Float64))
    cs = pop[T].size
    fs = pop[Female].size
    ms = pop[Male].size
    fwv = (weight(pop[Female].fit[1:fs[1]]), weight(pop[Female].fit[fs[1]+1:end]))
    mwv = (weight(pop[Male].fit[1:ms[1]]), weight(pop[Male].fit[ms[1]+1:end]))
    fr = (fb[1], 1 - fb[2])
    mr = (mb[1], 1 - mb[2])
    for i = 1:length(ids)
        deme = i <= cs[1] ? 1 : 2
        ids[i,1] = rand() < fr[deme] ? sample(fwv[1]) : sample(fwv[2]) + fs[1]
        ids[i,2] = rand() < mr[deme]? sample(mwv[1]) : sample(mwv[2]) + ms[1]
    end
    nothing
end

function pickparents!(ids::(ParentIds, ParentIds), pop::Population, fb::(Float64, Float64), mb::(Float64, Float64))
    pickparents!(Female, ids[1], pop, fb, mb)
    pickparents!(Male, ids[2], pop, fb, mb)
    nothing
end
