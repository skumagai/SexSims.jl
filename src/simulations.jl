using SexSims
using JSON
using Distributions: Binomial
using DataFrames: DataFrame, writetable
using DataArrays: DataArray, NA

# functions to parse input parameters from JSON and convert to appropriate format
function learningrate(ps, from, to)
    for elem in ps
        elem["from"] == from && elem["to"] == to && return convert(Vector{Float64}, elem["rate"])
    end
    Float64[0., 0.]
end

trait2fitness(ps) = [convert(Vector{Float64}, ps[1]) convert(Vector{Float64}, ps[2])]'

# core function for reproduction
function getchild(t, orgs, sex, deme, nf, nm, frate, mrate, sel, mut, rec)
    while true
        mom = rand(1:nf)
        pop = rand(1:nm)
        trait = SexSims.learn(deme, mom, pop, frate, mrate)
        if sel[deme, trait] < rand()
            return sex(t, orgs[Female].data[mom], orgs[Male].data[pop], mut, rec), trait
        end
    end
end

function setindividual!(orgs, idx, org, trait)
    orgs.data[idx] = org
    orgs.trait[idx] = trait
    nothing
end

function nmigrants!(mig, n, v)
    for i = 1:2
        mig[i] = rand(Binomial(n[i], v[i]))
    end
end

function nbefore!(nb, n, m)
    for i = 1:2
        nb[i] = n[i] + m[i] - m[3 - i]
    end
end

function simulation(config)
    params = JSON.parsefile(config)
    srand(1)
    tmax = 10
    # Population sizes of adults participating in reproduction.  These numbers are per-deme.
    nf = convert(Vector{Int}, params["population size"]["female"])
    nm = convert(Vector{Int}, params["population size"]["male"])
    # Proportion of reproducing migrants out of all reproducing adults.
    # This is sex- and deme-specific.
    fb = convert(Vector{Float64}, params["migration"]["female"]["fraction"])
    mb = convert(Vector{Float64}, params["migration"]["male"]["fraction"])
    # cost of migration.  This is sex- and deme-specific
    ft = convert(Vector{Float64}, params["migration"]["female"]["cost"])
    mt = convert(Vector{Float64}, params["migration"]["male"]["cost"])
    # Numbers of migrants participating in reproduction.  This quantity takes cost into account.
    # This is deme-specific and sex-specific in parents.
    fv = convert(Vector{Float64}, fb .* ft ./ (fb .* ft + 1 - fb))
    mv = convert(Vector{Float64}, mb .* mt ./ (mb .* mt + 1 - mb))
    # Probability of mutation per locus per generation.  This is common across loci.
    mut = convert(Float64, params["mutation probability"])
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

    # number of migrants
    fmig = Array(Int, 2)
    mmig = Array(Int, 2)
    # deme-specific number of offspring
    nfo = Array(Int, 2)
    nmo = Array(Int, 2)

    nmigrants!(fmig, nf, fv)
    nmigrants!(mmig, nm, mv)
    nbefore!(nfo, nf, fmig)
    nbefore!(nmo, nm, mmig)

    parents = Population(nfo, nmo)
    children = Population(nfo, nmo)
    rec = GeneStateRecorder(10 * (sum(nfo) + sum(nmo)))

    for t = 1:tmax
        parents, children = children, parents

        # mating
        fcidx = 1
        mcidx = 1
        for deme = 1:2
            for i = 1:nfo[deme]
                setindividual!(children[Female], fcidx, getchild(t, parents, Female, deme, nf[deme], nm[deme], f2fl[deme], m2fl[deme], ffit, mut, rec)...)
                fcidx += 1
            end
            for i = 1:nmo[deme]
                setindividual!(children[Male], mcidx, getchild(t, parents, Male, deme, nf[deme], nm[deme], f2ml[deme], m2ml[deme], mfit, mut, rec)...)
                mcidx += 1
            end
        end
        migrate!(t, children[Female], nfo, fmig, rec)
        migrate!(t, children[Male], nmo, mmig, rec)
        # decide number of local and migrating individuals in next-generation.
        nmigrants!(fmig, nf, fv)
        nmigrants!(mmig, nm, mv)
        nbefore!(nfo, nf, fmig)
        nbefore!(nmo, nm, mmig)
    end
    rec, children
end

function summarize(dir, rec, pop)
    getmigrations(joinpath(dir, "migrations.tsv"), rec, pop)
    getdistance(joinpath(dir, "dist.tsv"), rec, pop)
end

function listuniquegenes(genes)
    dict = Dict{SexSims.StateIndex, SexSims.Gene}()
    for gene in genes
        if !haskey(dict, SexSims.id(gene))
            dict[SexSims.id(gene)] = gene
        end
    end
    [elem for elem in values(dict)]
end

function recorddistance!(idx, rec, chrstr, chr, gene1, gene2, dist, alleles)
    for i = 1:length(alleles)
        for j = (i+1):length(alleles)
            chr[idx] = chrstr
            gene1[idx] = convert(Int, SexSims.id(alleles[i]))
            gene2[idx] = convert(Int, SexSims.id(alleles[j]))
            d = SexSims.distance(rec, alleles[i], alleles[j], SexSims.MUTATION)
            dist[idx] = isnull(d) ? NA : get(d)
            idx += 1
        end
    end
    idx
end

function getdistance(path, rec, pop)
    agene = vcat([[gene for gene in ind.auto] for ind in pop[SexSims.Female].data]...,
                 [[gene for gene in ind.auto] for ind in pop[SexSims.Male].data]...)
    xgene = vcat([[gene for gene in ind.auto] for ind in pop[SexSims.Female].data]...,
                 [ind.x for ind in pop[SexSims.Male].data])
    ygene = [ind.y for ind in pop[SexSims.Male].data]
    mitogene = [ind.mito for ind in pop[SexSims.Female].data]
    a = listuniquegenes(agene)
    x = listuniquegenes(xgene)
    y = listuniquegenes(ygene)
    mito = listuniquegenes(mitogene)
    nelem = binomial(length(a), 2) + binomial(length(x), 2) + binomial(length(y), 2) + binomial(length(mito), 2)
    chr = Array(ASCIIString, nelem)
    gene1 = Array(Int, nelem)
    gene2 = Array(Int, nelem)
    dist = DataArray(Int, nelem)
    idx = 1
    idx = recorddistance!(idx, rec, "autosome", chr, gene1, gene2, dist, a)
    idx = recorddistance!(idx, rec, "x", chr, gene1, gene2, dist, x)
    idx = recorddistance!(idx, rec, "y", chr, gene1, gene2, dist, y)
    idx = recorddistance!(idx, rec, "mito", chr, gene1, gene2, dist, mito)
    data = DataFrame(chromosome = chr, allele_id_1 = gene1, allele_id_2 = gene2, distance = dist)
    writetable(path, data, separator = '\t')
end

function getmigrations(path,rec, pop)
    len = 5 * length(pop[SexSims.Female].data) + 4 * length(pop[SexSims.Male].data)
    ids = Array(Int, len)
    demes = Array(Int, len)
    sexes = Array(ASCIIString, len)
    chrs = Array(ASCIIString, len)
    alleles = Array(Int, len)
    nmigs = Array(Int, len)

    idx = 1
    i = 1
    d = 1
    loc = 1
    s = pop[SexSims.Female].size[1]
    for org in pop[SexSims.Female].data
        loc = d <= s ? 1 : 2
        for (c, chrstr) in zip((:auto, :x), ("autosome", "x"))
            for l in getfield(org, c)
                ids[idx] = i
                sexes[idx] = "female"
                chrs[idx] = chrstr
                alleles[idx] = convert(Int, SexSims.id(l))
                nmigs[idx] = SexSims.countalong(rec, l, SexSims.MIGRATION)
                demes[idx] = loc
                idx += 1
            end
        end
        ids[idx] = i
        sexes[idx] = "female"
        chrs[idx] = "mitochondrion"
        alleles[idx] = convert(Int, SexSims.id(org.mito))
        nmigs[idx] = SexSims.countalong(rec, org.mito, SexSims.MIGRATION)
        demes[idx] = loc
        idx += 1
        i += 1
        d += 1
    end
    d = 1
    s = pop[SexSims.Male].size[1]
    for org in pop[SexSims.Male].data
        loc = d <= s ? 1 : 2
        for a in org.auto
            ids[idx] = i
            sexes[idx] = "male"
            chrs[idx] = "autosome"
            alleles[idx] = convert(Int, SexSims.id(a))
            nmigs[idx] = SexSims.countalong(rec, a, SexSims.MIGRATION)
            demes[idx] = s
            idx += 1
        end
        for (c, chrstr) in zip((:x, :y), ("x", "y"))
            ids[idx] = i
            sexes[idx] = "male"
            chrs[idx] = chrstr
            alleles[idx] = convert(Int, SexSims.id(getfield(org, c)))
            nmigs[idx] = SexSims.countalong(rec, getfield(org, c), SexSims.MIGRATION)
            demes[idx] = s
            idx += 1
        end
        i += 1
        d += 1
    end

    data = DataFrame(organism_id = ids, sex = sexes, deme = demes ,chromosome = chrs, allele_id = alleles , number_of_migrations = nmigs)
    writetable(path, data, separator = '\t')
end

function getresultdir(config)
    i = 1
    while isdir("$config-$i")
        i += 1
    end
    dir = "$config-$i"
    mkdir(dir)
    dir
end

function record(path, str)
    f = open(path, "a")
    write(f, str)
    close(f)
end

function main()
    length(ARGS) != 1 && error("Usage: julia simulation.jl inputfile")
    config = ARGS[1]
    isfile(config) || error("Not a file: $config")
    resultdir = getresultdir(config)
    logf = joinpath(resultdir, "simulation.log")
    seed = rand(1:typemax(UInt32))
    record(logf, "$(strftime(time())) Simulation launched\n")
    record(logf, "$(strftime(time())) Results stored at: $resultdir\n")
    record(logf, "$(strftime(time())) RNG Seed: $(rand(1:seed))\n")
    srand(seed)
    rec, pop = simulation(config)
    record(logf, "$(strftime(time())) Simulation finished\n")
    record(logf, "$(strftime(time())) Saving results\n")
    summarize(resultdir, rec, pop)
    record(logf, "$(strftime(time())) Done\n")
end

main()
