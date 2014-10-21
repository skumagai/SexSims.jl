using StatsBase: WeightVec

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

    ft = vcat([fill(convert(DemeIndex, i), j) for (i, j) = enumerate(nf)]...)
    mt = vcat([fill(convert(DemeIndex, i), j) for (i, j) = enumerate(nm)]...)

    f = Organisms{Female}(nf, ffit, ft, [Female() for _ = 1:sum(nf)])
    m = Organisms{Male}(nm, mfit, mt, [Male() for _ = 1:sum(nm)])

    Population(f, m)
end

Base.getindex(p::Population, ::Type{Female}) = p.f
Base.getindex(p::Population, ::Type{Male}) = p.m

immutable ParentIds
    data::Matrix{Int}
end
ParentIds(n::(Int, Int)) = ParentIds(Array(Int, sum(n), 2))
ParentIds(n::Int) = ParentIds(Array(Int, n, 2))

Base.getindex(p::ParentIds, i) = Base.getindex(p.data,  i, 1:2)
Base.getindex(p::ParentIds, i, j) = Base.getindex(p.data, i, j)
Base.setindex!(p::ParentIds, i, j, k) = Base.setindex!(p.data, i, j, k)

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

function fitness!{T<:Organism}(p::Organisms{T}, r1::Float64, r2::Float64)
    for i = 1:length(p.data)
        deme = i <= p.size[1] ? 1 : 2
        r = deme == 1 ? r1 : r2
        p.fit[i] = fitness(deme, p.trait[i], 1.0, r)
    end
    nothing
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

function pickparents!(col, ids::ParentIds, w1::WeightVec, w2::WeightVec, r::(Float64, Float64))
    offset = 0
    for (deme, w) = enumerate((w1, w2))
        for i = (1:length(w)) + offset
            ids[i, col] = pickparent(w1, w2, r[deme])
        end
        offset += length(w)
    end
    nothing
end

function pickparents!(size::(Int, Int), ids::ParentIds, pop::Population, fb::(Float64, Float64), mb::(Float64, Float64))
    fwv = (WeightVec(pop[Female].fit[1:size[1]]), WeightVec(pop[Female].fit[size[1]+1:end]))
    pickparents!(1, ids, fwv..., (fb[1], 1 - fb[2]))

    mwv = (WeightVec(pop[Male].fit[1:size[1]]), WeightVec(pop[Male].fit[size[1]+1:end]))
    pickparents!(2, ids, mwv..., (mb[1], 1 - mb[2]))
    nothing
end
