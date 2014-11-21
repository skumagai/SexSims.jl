using StatsBase: WeightVec, sample

immutable Organism{T<:(Chromosome...)}
    chrs::T
end

makeancestor(chrlist) = Organism(tuple([makeancestor(chr) for chr in chrlist]...))

function Organism(t, mom, pop, mut, rec, sex::Type{Sex})
    dchr = [fuse(meiosis(t, mchr, mut, rec, Female, sex),
                 meiosis(t, pchr, mut, rec, Male, sex))
            for (mchr, pchr) in zip(mom.chrs, pop.chrs)]

    Orgnasims(tuple(dchr...))
end

Base.length(o::Organism) = length(o.chrs)

function migrate!{S<:Sex}(t, p, rec, sex::Type{S})
    Organism(tuple([migrate!(t, c, rec, sex) for c in p.chrs]...))
end
