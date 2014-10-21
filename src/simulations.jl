using StatsBase: sample

function simulation()
    srand(1234)
    tmax = 10
    nf = (100, 200)
    nm = (200, 100)
    fb = (0.2, 0.3)
    mb = (0.3, 0.2)
    mut = 0.02
    fsel = (0.1, 0.2)
    msel = (0.2, 0.2)

    parents = Population(nf, nm)
    children = Population(nf, nm)
    rec = GeneStateRecorder(4 * (sum(nf) + sum(nm)))

    # Index of parents
    pids = (ParentIds{Female}(nf), ParentIds{Male}(nm))

    for t = 1:tmax
        parents, offspring = offspring, parents
        # mating
        pickparents!(pids, parents, fb, mb)
        reproduce!(t, children, parents, pids, mut, rec)
        fitness!(children, fsel, msel)
        # learning
        learn!(children, parents, pids, ftof, mtof, ftom, mtom)
    end
    rec
end
