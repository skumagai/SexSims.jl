using SexSims
using JSON

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

    # Index of parents
    fpids = ParentIds(nf)
    mpids = ParentIds(nm)

    for t = 1:tmax
        parents, children = children, parents
        fitness!(parents[Female], ffit...)
        fitness!(parents[Male], mfit...)
        # mating
        # Select aprents of female offspring
        pickparents!(nf, fpids, parents, fv, mv)
        # Select parents of male offspring
        pickparents!(nm, mpids, parents, fv, mv)
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
