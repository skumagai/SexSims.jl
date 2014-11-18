using SexSims
using JSON
using Distributions: Binomial
using DataFrames: DataFrame, writetable

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

function summarize(config, rec, pop)
    i = 1
    while isfile("$config-$i")
        i += 1
    end
    outfile = "$config-$i"
    mkdir(outfile)
    getmigrations(joinpath(outfile, "migrations.tsv"))
    getdistance(joinpath(outfile, "dist.tsv"))
end

function getdistance(path, rec, pop)
    a = unique(vcat(unique([SexSims.id(gene) for ind in pop[SexSims.Female] for gene in ind.auto]),
                    unique([SexSims.id(gene) for ind in pop[SexSims.Male] for gene in ind.auto])))
    x = unique(vcat(unique([SexSims.id(gene) for ind in pop[SexSims.Female] for gene in ind.x]),
                    unique([SexSims.id(ind.x) for ind in pop[SexSims.Male]])))
    y = unique([SexSims.id(ind.y) for ind in pop[SexSims.Male]])
    mito = unique([SexSims.id(ind.mito) for ind.pop[SexSims.Female]])
    sort!(a)
    sort!(x)
    sort!(y)
    sort!(mito)
    nelem = length(a) * (length(a) - 1) + length(x) * (length(x) - 1) + length(y) * (length(y) - 1) + length(mito) * (length(mito) - 1)
    chr = Array(ASCIIString, nelem)
    gene1 = Array(Int, nelem)
    gene2 = Array(Int, nelem)
    dist = DataArray(Int, nelem)
    idx = 1
    for i = 1:length(a)
        for j = (i+1):length(a)
            chr[idx] = "autosome"
            gene1[idx] = SexSims.id(a[i])
            gene2[idx] = SexSims.id(a[j])
            d = SexSims.distance(rec, mito[i], mito[j], "mutation")
            dist[idx] = isnull(d) ? NA : get(d)
            idx += 1
        end
    end
    for i = 1:length(x)
        for j = (i+1):length(x)
            chr[idx] = "autosome"
            gene1[idx] = SexSims.id(x[i])
            gene2[idx] = SexSims.id(x[j])
            d = SexSims.distance(rec, mito[i], mito[j], "mutation")
            dist[idx] = isnull(d) ? NA : get(d)
            idx += 1
        end
    end
    for i = 1:length(y)
        for j = (i+1):length(y)
            chr[idx] = "autosome"
            gene1[idx] = SexSims.id(y[i])
            gene2[idx] = SexSims.id(y[j])
            d = SexSims.distance(rec, mito[i], mito[j], "mutation")
            dist[idx] = isnull(d) ? NA : get(d)
            idx += 1
        end
    end
    for i = 1:length(mito)
        for j = (i+1):length(mito)
            chr[idx] = "autosome"
            gene1[idx] = SexSims.id(mito[i])
            gene2[idx] = SexSims.id(mito[j])
            d = SexSims.distance(rec, mito[i], mito[j], "mutation")
            dist[idx] = isnull(d) ? NA : get(d)
            idx += 1
        end
    end
    data = DataFrame(chr = chr, gene1 = gene1, gene2 = gene2, distance = dist)
    writetable(path, data, separator = '\t')
end

function getmigrations(path)

    println("deme\tsex\torganism\tchromosome\tallele.index\tmigration.index")
    orgs = pop[Female]
    b = orgs.size[1]
    idx = 1
    for i = 1:length(orgs.data)
        deme = i <= b ? 1 : 2
        for j = 1:2
            println("$deme\tfemale\t$idx\tautosome\t$(SexSims.id(orgs.data[i].auto[j]))\t$(SexSims.mig(orgs.data[i].auto[j]))")
        end
        for j = 1:2
            println("$deme\tfemale\t$idx\tx\t$(SexSims.id(orgs.data[i].x[j]))\t$(SexSims.mig(orgs.data[i].x[j]))")
        end
        println("$deme\tfemale\t$idx\tmito\t$(SexSims.id(orgs.data[i].mito))\t$(SexSims.mig(orgs.data[i].mito))")
        idx += 1
    end
    orgs = pop[Male]
    b = orgs.size[1]
    for i = 1:length(orgs.data)
        deme = i <= b ? 1 : 2
        for j = 1:2
            println("$deme\tmale\t$idx\tautosome\t$(SexSims.id(orgs.data[i].auto[j]))\t$(SexSims.mig(orgs.data[i].auto[j]))")
        end
        println("$deme\tmale\t$idx\tx\t$(SexSims.id(orgs.data[i].x))\t$(SexSims.mig(orgs.data[i].x))")
        println("$deme\tmale\t$idx\ty\t$(SexSims.id(orgs.data[i].y))\t$(SexSims.mig(orgs.data[i].y))")
        idx += 1
    end
end

function main()
    length(ARGS) != 1 && error("Usage: julia simulation.jl inputfile")
    config = ARGS[1]
    isfile(config) || error("Not a file: $config")
    rec, pop = simulation(config)
    summarize(config, rec, pop)
end

main()
