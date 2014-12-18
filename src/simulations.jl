using SexSims: migrate!, getresultdir, openlog, writelog, closelog, nmigrants!, nbefore!, Parameters, Organism, Autosome, XChromosome, YChromosome, Mitochondrion, Female, Male, makeancestor, GeneStateRecorder, learn, savedistance, savemigrations, getallgenes, samplemigintervals
using StatsBase: sample

@inline indeme1(x) = 0x1 == x ? true : false
@inline indeme2(x) = !indeme1(x)

function simulate(config)
    p = Parameters(config)

    # number of migrants
    migs = (Array(Int, 2), Array(Int, 2))
    # deme-specific number of offspring
    no = (Array(Int, 2), Array(Int, 2))

    for i = 1:2
        nmigrants!(migs[i], p.n[i], p.v[i])
        nbefore!(no[i], p.n[i], migs[i])
    end

    ps = (Array(Organism{(Autosome{1}, XChromosome{1}, YChromosome{1}, Mitochondrion{1})}, sum(p.n[1])),
          Array(Organism{(Autosome{1}, XChromosome{1}, YChromosome{1}, Mitochondrion{1})}, sum(p.n[2])))
    cs = ([makeancestor((Autosome{1}, XChromosome{1}, YChromosome{1}, Mitochondrion{1}))
           for _ = 1:sum(p.n[1])],
          [makeancestor((Autosome{1}, XChromosome{1}, YChromosome{1}, Mitochondrion{1}))
           for _ = 1:sum(p.n[2])])

    ploc = (Array(UInt8, sum(p.n[1])), Array(UInt8, sum(p.n[2])))
    cloc = (vcat(ones(UInt8, p.n[1][1]), 2 * ones(UInt8, p.n[1][2])),
            vcat(ones(UInt8, p.n[2][1]), 2 * ones(UInt8, p.n[2][2])))
    ptrait = deepcopy(ploc)
    ctrait = deepcopy(cloc)

    rec = GeneStateRecorder(20 * (sum(p.n[1]) + sum(p.n[2])))

    indeme = (indeme1, indeme2)

    for t = 1:p.tmax
        t % 1000 == 0 && info("generation: $t")
        ps, cs = cs, ps
        ptrait, ctrait = ctrait, ptrait
        ploc, cloc = cloc, ploc

        # mating
        for (c, loc, trait, sex, j) in zip(cs, cloc, ctrait, (Female, Male), 1:2)
            # decide number of local and migrating individuals in next-generation.
            nmigrants!(migs[j], p.n[j], p.v[j])
            nbefore!(no[j], p.n[j], migs[j])

            idx = 1
            for deme = 0x1:0x2
                mpos = find(indeme[deme], ploc[1])
                ppos = find(indeme[deme], ploc[2])
                miglist = rand(1:p.n[j][deme], migs[j][deme])
                for i = 1:no[j][deme]
                    mid = sample(mpos)
                    pid = sample(ppos)
                    tval = learn(deme, ptrait[1][mid], ptrait[2][pid], p.l[1][j][deme], p.l[2][j][deme])
                    while rand() > p.fit[j][deme, tval]
                        mid = sample(mpos)
                        pid = sample(ppos)
                        tval = learn(deme, ptrait[1][mid], ptrait[2][pid], p.l[1][j][deme], p.l[2][j][deme])
                    end
                    trait[idx] = tval
                    if i in miglist
                        loc[idx] = 0x3 - deme
                        o = Organism(t, ps[1][mid], ps[2][pid], p.mut, rec, sex)
                        c[idx] = migrate!(t, o, rec, sex)
                    else
                        loc[idx] = deme
                        c[idx] = Organism(t, ps[1][mid], ps[2][pid], p.mut, rec, sex)
                    end
                    idx += 1
                end
            end
        end
    end
    rec, cs, cloc, ctrait
end

function summarize(dir, rec, pop, demes, traits)
    df = getallgenes(pop, demes)
    savemigrations(joinpath(dir, "migcounts.tsv"), rec, df)
    savedistance(joinpath(dir, "dist.tsv"), rec, df)
    samplemigintervals(joinpath(dir, "migintervals.tsv"), rec, df)
    savefst(joinpath(dir, "fst.tsv"), rec, df)
end


function main()
    length(ARGS) != 1 && error("Usage: julia simulation.jl inputfile")
    config = ARGS[1]
    isfile(config) || error("Not a file: $config")
    resultdir = getresultdir(config)
    logfile = joinpath(resultdir, "simulation.log")
    logh = openlog(logfile)
    seed = rand(1:typemax(UInt32))
    writelog(logh, "$(strftime(time())) Simulation launched\n")
    @osx_only writelog(logh, readall(`sysctl kern.version kern.hostname hw.model hw.ncpu hw.physmem hw.cpufamily hw.cputype hw.cpusubtype`))
    writelog(logh, "$(strftime(time())) Results stored in: $resultdir\n")
    writelog(logh, "$(strftime(time())) RNG Seed: $(rand(1:seed))\n")
    srand(seed)
    rec, pop, demes, traits = simulate(config)
    writelog(logh, "$(strftime(time())) Simulation finished\n")
    writelog(logh, "$(strftime(time())) Saving results\n")
    summarize(resultdir, rec, pop, demes, traits)
    writelog(logh, "$(strftime(time())) Done\n")
    closelog(logh)
end

main()
