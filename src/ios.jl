using JSON

# Input functions
immutable Parameters
    # 1st element is for female-parent, and 2nd element is for male-parent.
    # Locus population sizes right before mating takes place.
    n::NTuple{2, Vector{Int}}
    # Proportion of reproducing migrants out of all reproducing adults corrected by
    # cost of migration.
    v::NTuple{2, Vector{Float64}}
    # Mutation rate
    mut::Float64
    # Fitness of non-local organisms (deme-by-trait)
    fit::NTuple{2, Matrix{Float64}}
    # Learning rates
    l::NTuple{2, NTuple{2, Vector{Float64}}}
end

function readinput(config)
    params = JSON.parsefile(config)
    nf = convert(Vector{Int}, params["population size"]["female"])
    nm = convert(Vector{Int}, params["population size"]["male"])
    fb = convert(Vector{Float64}, params["migration"]["female"]["fraction"])
    mb = convert(Vector{Float64}, params["migration"]["male"]["fraction"])
    ft = convert(Vector{Float64}, params["migration"]["female"]["cost"])
    mt = convert(Vector{Float64}, params["migration"]["male"]["cost"])
    fv = convert(Vector{Float64}, fb .* ft ./ (fb .* ft + 1 - fb))
    mv = convert(Vector{Float64}, mb .* mt ./ (mb .* mt + 1 - mb))
    mut = convert(Float64, params["mutation probability"])
    ffit = trait2fitness(params["trait"]["female"])
    mfit = trait2fitness(params["trait"]["male"])
    f2fl = learningrate(params["learning"], "female", "female")
    m2fl = learningrate(params["learning"], "male", "female")
    f2ml = learningrate(params["learning"], "female", "male")
    m2ml = learningrate(params["learning"], "male", "male")
    Parameters((nf, nm), (fv, mv), mut, (ffit, mfit), ((f2fl, f2ml), (m2fl, m2ml)))
end

function learningrate(ps, from, to)
    for elem in ps
        elem["from"] == from && elem["to"] == to && return convert(Vector{Float64}, elem["rate"])
    end
    Float64[0., 0.]
end

trait2fitness(ps) = [convert(Vector{Float64}, ps[1]) convert(Vector{Float64}, ps[2])]'

# Output functions
function getresultdir(config)
    i = 1
    while isdir("$config-$i")
        i += 1
    end
    dir = "$config-$i"
    mkdir(dir)
    dir
end

openlog(logfile) = open(logfile, "a")
closelog(log) = close(log)
writelog(log, str) = write(log, str)

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

function savedistance(path, rec, pop)
    n = [length(pop[1]), length(pop[2])]
    #
    nloci = 0
    ngenes = 0
    chrlist = typeof(pop[1][1])

    for chr in chrlist
        if chr <: Autosome
            ngenes += 2 * chr.parameters[1] * sum(n)
            nloci += 1
        elseif chr <: XYChromosome
            ngenes += 2 * chr.parameters[1] * sum(n)
            nloci += 2
        else

            ngenes += chr.parameters[1] * n[1]
        end
    end

    genes = Array(Gene, ngenes)
    loci = Array(Int, ngenes)
    inds = Array(Int, ngenes)

    colidx = 1
    indidx = 1
    for org in pop[1]
        locusidx = 1
        for c in p.chrs, gs in getgenotype(c)
            if isa(c, Autosome)
                for g in gs
                    genes[colidx] = g
                    chrnames[colidx] = "autosome"
                    loci

            elseif isa(c, XYChromosome)

            else
            end
            colidx += 1
            locusidx += 1
        end
        indidx += 1
    end


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
    data = DataFrame(chromosome = chr, locus_id = loci, allele_id_1 = gene1, allele_id_2 = gene2, distance = dist)
    writetable(path, data, separator = '\t')
end

function savemigrations(path,rec, pop)
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
