using StatsBase: sample!
using SexSims
using JSON

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

    # Index of parents
    fpids = ParentIds(nfo)
    mpids = ParentIds(nmo)

    fmigs = [Array(Int, v) for v in fv]
    mmigs = [Array(Int, v) for v in mv]

    for t = 1:tmax
        parents, children = children, parents
        # migration
        # select migrants
        samples!(1:fv[1], fmigs[1], replace = false)
        samples!(nfo[1]+(1:fv[2]), fmigs[2], replace = false)
        samples!(1:mv[1], mmigs[1], replace = false)
        samples!(nmo[1]+(1:mv[2]), mmigs[2], replace = false)


        # mating
        # Select aprents of female offspring
        pickparents!(nfo, fpids, parents)
        # Select parents of male offspring
        pickparents!(nmo, mpids, parents)

        fitness!(parents[Female], ffit...)
        fitness!(parents[Male], mfit...)
        # Generating offspring
        reproduce!(t, children, parents, (fpids, mpids), mut, rec)
        # learning
        learn!(children[Female], parents, fpids, f2fl, m2fl)
        learn!(children[Male], parents, mpids, f2ml, m2ml)
    end
    rec, children
end

function summarize(rec, pop)
end

function main()
    length(ARGS) != 1 && error("Usage: julia simulation.jl inputfile")
    simulation()
end

main()
