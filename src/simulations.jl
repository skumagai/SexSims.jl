using StatsBase: sample!
using SexSims
using JSON
using StatsBase: sample
using Distributions: Binomial

function learningrate(ps, from, to)
    for elem in ps
        elem["from"] == from && elem["to"] == to && return convert(Array{Float64,1}, elem["rate"])
    end
    Float64[0., 0.]
end

trait2fitness(ps) = Float64[ps[1][2], ps[2][1]]

function simulation()
    params = JSON.parsefile(ARGS[1])
    srand(1)
    tmax = 10
    # Population sizes of adults participating in reproduction.  These numbers are per-deme.
    nf = convert(Array{Int,1}, params["population size"]["female"])
    nm = convert(Array{Int,1}, params["population size"]["male"])
    # Proportion of reproducing migrants out of all reproducing adults.
    # This is sex- and deme-specific.
    fb = convert(Array{Float64,1}, params["migration"]["female"]["fraction"])
    mb = convert(Array{Float64,1}, params["migration"]["male"]["fraction"])
    # cost of migration.  This is sex- and deme-specific
    ft = convert(Array{Float64,1}, params["migration"]["female"]["cost"])
    mt = convert(Array{Float64,1}, params["migration"]["male"]["cost"])
    # Numbers of migrants participating in reproduction.  This quantity takes cost into account.
    # This is deme-specific and sex-specific in parents.
    fv = Int[round(val) for val in nf .* fb .* ft ./ (fb .* ft + 1 - fb)]
    mv = Int[round(val) for val in nm .* mb .* mt ./ (mb .* mt + 1 - mb)]
    # number of offspring before migration
    nfo = nf - reverse(fv) + fv
    nmo = nm - reverse(mv) + mv
    # Probability of mutation per locus per generation.  This is common across loci.
    mut = params["mutation probability"]
    # Fitness of individuals carrying non-local trait (relative to local individuals)
    # This is deme specific (in deme 1, in deme 2)
    ffit = trait2fitness(params["trait"]["female"])
    mfit = trait2fitness(params["trait"]["male"])
    # Probabilities of learning parental trait.
    # This is sex-specific in both parents and offspring as well as deme-specific.
    f2fl = learningrate(params["learning"], "female", "female")
    m2fl = learningrate(params["learning"], "male", "female")
    f2ml = learningrate(params["learning"], "female", "male")
    m2ml = learningrate(params["learning"], "male", "male")

    parents = Population(nfo, nmo)
    children = Population(nfo, nmo)
    rec = GeneStateRecorder(4 * (sum(nfo) + sum(nmo)))

    nfo = Array(Int, 2)
    nmo = Array(Int, 2)
    mfo = Array(Int, 2)
    mmo = Array(Int, 2)
    for t = 1:tmax
        parents, children = children, parents
        # decide number of local and migrating offspring.
        for i = 1:2
            mfo[i] = rand(Binomial(nf[i], fv[i]))
            mmo[i] = rand(Binomial(nm[i], mv[i]))
        end
        for i = 1:2
            nfo[i] = nf[i] - mfo[i] + mfo[3 - i]
            nmo[i] = nm[i] - mmo[i] + mmo[3 - i]
        end

        # mating
        fc = 1
        mc = 1
        for deme = 1:2
            for i = 1:nfo[deme]
                while true
                    mom = rand(1:nf[deme])
                    pop = rand(1:mn[deme])
                    trait = learn(deme, mom, pop, f2fl, m2fl)
                    if sel[deme, trait] < rand()
                        children[Female][fc] = Female(mom, pop, trait)
                        break
                    end
                end
                fc += 1
            end
            for i = 1:nmo[deme]
                while true
                    mom = rand(1:nf[deme])
                    pop = rand(1:nm[deme])
                    trait = learn(deme, mom, pop, f2ml, m2ml)
                    if sel[deme, trait] < rand()
                        children[Male][mc] = Male(mom, pop, trait)
                        break
                    end
                end
                mc += 1
            end
        end
        migrate!(t, femeles, nfo, mfo, rec)
        migrate!(t, males, nmo, mmo, rec)
    end
    rec, children
end

function migrate!(t, orgs, n, migs, rec)
    tar, src = getmigrants(no, migs)
    b = n[1]
    for i = 1:length(tar)
        if tar[i] <= b < src[i] || src[i] <= b < tar[i]
            orgs[src[i]] = migrate!(t, orgs[src[i]], rec)
        end
    end
    orgs[tar] = orgs[src]
end



function summarize(rec, pop)
end

function main()
    length(ARGS) != 1 && error("Usage: julia simulation.jl inputfile")
    simulation()
end

main()
