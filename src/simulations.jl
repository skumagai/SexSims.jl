using SexSims: readinput, migrate!, getresultdir, openlog, writelog, closelog, nmigrants!, nbefore!, Organism, Autosome, XYChromosome, Mitochondrion, Female, Male, makeancestor, GeneStateRecorder, learn
using JSON
using Distributions: Binomial
using DataFrames: DataFrame, writetable
using DataArrays: DataArray, NA

@inline indeme1(x) = 0x1 == x ? true : false
@inline indeme2(x) = !deme1(x)

function simulate(config)
    p = readinput(config)
    srand(1)
    tmax = 10

    # number of migrants
    fmig = Array(Int, 2)
    mmig = Array(Int, 2)
    # deme-specific number of offspring
    nfo = Array(Int, 2)
    nmo = Array(Int, 2)

    nmigrants!(fmig, p.n[1], p.v[1])
    nmigrants!(mmig, p.n[2], p.v[2])
    nbefore!(nfo, p.n[1], fmig)
    nbefore!(nmo, p.n[2], mmig)

    ps = (Array(Organism{(Autosome{1}, XYChromosome{1}, Mitochondrion{1})}, sum(p.n[1])),
          Array(Organism{(Autosome{1}, XYChromosome{1}, Mitochondrion{1})}, sum(p.n[2])))
    cs = ([makeancestor((Autosome{1}, XYChromosome{1}, Mitochondrion{1}))
           for _ = 1:sum(p.n[1])],
          [makeancestor((Autosome{1}, XYChromosome{1}, Mitochondrion{1}))
           for _ = 1:sum(p.n[2])])

    ploc = (Array(UInt8, sum(p.n[1])), Array(UInt8, sum(p.n[2])))
    cloc = (vcat(ones(UInt8, p.n[1][1]), 2 * ones(UInt8, p.n[1][2])),
            vcat(ones(UInt8, p.n[2][1]), 2 * ones(UInt8, p.n[2][2])))
    ptrait = deepcopy(ploc)
    ctrait = deepcopy(cloc)

    rec = GeneStateRecorder(10 * (sum(p.n[1]) + sum(p.n[2])))

    indeme = (indeme1, indeme2)

    for t = 1:tmax
        println(t)
        ps, cs = cs, ps
        ptrait, ctrait = ctrait, ptrait
        ploc, cloc = cloc, ploc

        # mating
        for (c, loc, trait, n, m, sex, j) in zip(cs, cloc, ctrait, (nfo, nmo), (fmig, mmig), (Female, Male), 1:2)
            idx = 1
            for deme = 0x1:0x2
                mpos = find(indeme[deme], ploc[1])
                ppos = find(indeme[deme], ploc[2])
                migs = rand(1:n[deme], m[deme])
                for i = 1:n[deme]
                    mid = rand(mpos)
                    pid = rand(ppos)
                    tval = learn(deme, ptrait[1][mid], ptrait[2][pid], p.l[1][j][deme], p.l[2][j][deme])
                    while p.fit[j][deme, t] > rand()
                        mid = rand(mpos)
                        pid = rand(ppos)
                        tval = learn(deme, ptrait[1][mid], ptrait[2][pid], p.l[1][j][deme], p.l[2][j][deme])
                    end
                    trait[idx] = tval
                    if idx in migs
                        loc[idx] = 0x3 - deme
                        o = Organism(t, ps[1][mid], ps[2][pid], p.mut, rec, sex)
                        c[idx] = migrate!(t, o, rec, sex)
                    else
                        loc[idx] = deme
                        c[idx] = Organism(t, ps[1][mid], ps[2][pid], p.mut, rec, sex)
                    end
                    didx += 1
                end
            end
        end
        # decide number of local and migrating individuals in next-generation.
        nmigrants!(fmig, p.n[1], p.v[1])
        nmigrants!(mmig, p.n[2], p.v[2])
        nbefore!(nfo, p.n[1], fmig)
        nbefore!(nmo, p.n[2], mmig)
    end
    rec, cs
end

function summarize(dir, rec, pop)
    getmigrations(joinpath(dir, "migrations.tsv"), rec, pop)
    getdistance(joinpath(dir, "dist.tsv"), rec, pop)
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
    rec, pop = simulate(config)
    writelog(logh, "$(strftime(time())) Simulation finished\n")
    writelog(logh, "$(strftime(time())) Saving results\n")
#    summarize(resultdir, rec, pop)
    writelog(logh, "$(strftime(time())) Done\n")
    closelog(logh)
end

main()
