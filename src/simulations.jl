using SexSims
using JSON
using StatsBase: sample
using Distributions: Binomial

function simulation()
    params = JSON.parsefile(ARGS[1])
    println(params)
    srand(1)
    tmax = 10
    nf = (100, 200)
    nm = (200, 100)
    # Proportion of reproducing migrants out of all reproducing adults.
    # This is sex- and deme-specific.
    fb = [0.1, 0.3]
    mb = [0.2, 0.1]
    # cost of migration.  This is sex- and deme-specific
    ft = [0.2, 0.3]
    mt = [0.4, 0.1]
    # Actual proportion of reproducing migrants.  This quantity takes cost into account.
    # This is deme-specific and sex-specific in parents.
    fv = tuple(fb .* ft ./ (fb .* ft + 1 - fb)...)
    mv = tuple(mb .* mt ./ (mb .* mt + 1 - mb)...)
    # Probability of mutation per locus per generation.  This is common across loci.
    mut = 0.02
    # Fitness of individuals carrying non-local trait (relative to local individuals)
    # This is deme specific (in deme 1, in deme 2)
    ffit = (0.2, 0.3)
    mfit = (0.8, 0.9)
    # Probabilities of learning parental trait.
    # This is sex-specific in both parents and offspring as well as deme-specific.
    f2fl = (0.1, 0.2)
    m2fl = (0.2, 0.3)
    f2ml = (0.4, 0.2)
    m2ml = (0.2, 0.1)

    parents = Population(nf, nm)
    children = Population(nf, nm)
    rec = GeneStateRecorder(4 * (sum(nf) + sum(nm)))

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
