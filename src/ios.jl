using JSON
using DataArrays: DataArray, NA
using DataFrames: DataFrame, writetable

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
    # Number of generations to simulate
    tmax::UInt
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
    mut = float64(params["mutation probability"])
    ffit = trait2fitness(params["trait"]["female"])
    mfit = trait2fitness(params["trait"]["male"])
    f2fl = learningrate(params["learning"], "female", "female")
    m2fl = learningrate(params["learning"], "male", "female")
    f2ml = learningrate(params["learning"], "female", "male")
    m2ml = learningrate(params["learning"], "male", "male")
    t = float64(params["tmax"])
    Parameters((nf, nm), (fv, mv), mut, (ffit, mfit), ((f2fl, f2ml), (m2fl, m2ml)), t)
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

function getdistance!(df, idx, data, rec)
    chr, locus, genes = string(data[1]), data[2], data[3]
    for i = 1:length(genes)
        for j = (i+1):length(genes)
            df[idx, :chr_type] = chr
            df[idx, :locus_id] = locus
            df[idx, :gene_id_1] = genes[i]
            df[idx, :gene_id_2] = genes[j]
            d = distance(rec, genes[i], genes[j])
            df[idx, :distance] = isnull(d) ? NA : get(d)
            idx += 1
        end
    end
    idx
end

function nelems(pop)
    n = [length(pop[1]), length(pop[2])]
    ngenes = 0

    for chr in typeof(pop[1][1]).parameters[1]
        ntmp = chr.parameters[1]
        if chr <: Autosome
            ngenes += 2 * ntmp * sum(n)
        elseif chr <: XChromosome
            ngenes += ntmp * (2 * n[1] + n[2])
        elseif chr <: YChromosome
            ngenes += ntmp * n[2]
        elseif chr <: Mitochondrion
            ngenes += ntmp * n[1]
        end
    end
    ngenes
end

function getchrsymbol(chr)
    if isa(chr, Autosome)
        :autosome
    elseif isa(chr, XChromosome)
        :xchromosome
    elseif isa(chr, YChromosome)
        :ychromosome
    elseif isa(chr, Mitochondrion)
        :mitochondrion
    else
        :na
    end
end

function setgene!(df, colidx, indidx, sex, deme, chr, locus, gene)
    df[:org_id][colidx] = indidx
    df[:sex][colidx] = sex
    df[:deme][colidx] = deme
    df[:chr_type][colidx] = chr
    df[:locus_id][colidx] = locus
    df[:gene][colidx] = gene
    df
end

function getgenes!(pop, demes, df, colidx, indidx, ::Type{Female})
    for (ind, d) in zip(pop, demes)
        locidx = 1
        for chr in ind
            c = getchrsymbol(chr)
            if isa(chr, Autosome) || isa(chr, XChromosome)
                for i in 1:length(chr.loci1)
                    setgene!(df, colidx, indidx, :female, d, c, locidx, chr.loci1[i])
                    colidx += 1
                    setgene!(df, colidx, indidx, :female, d, c, locidx, chr.loci2[i])
                    colidx += 1
                    locidx += 1
                end
            elseif isa(chr, YChromosome)
                locidx += length(chr.loci1)
            elseif isa(chr, Mitochondrion)
                for i = 1:length(chr.loci1)
                    setgene!(df, colidx, indidx, :female, d, c, locidx, chr.loci1[i])
                    locidx += 1
                    colidx += 1
                end
            end
        end
        indidx += 1
    end
    colidx, indidx
end

function getgenes!(pop, demes, df, colidx, indidx, ::Type{Male})
    for (ind, d) in zip(pop, demes)
        locidx = 1
        for chr in ind
            c = getchrsymbol(chr)
            if isa(chr, Autosome)
                for i in 1:length(chr.loci1)
                    setgene!(df, colidx, indidx, :male, d, c, locidx, chr.loci1[i])
                    colidx += 1
                    setgene!(df, colidx, indidx, :male, d, c, locidx, chr.loci2[i])
                    colidx += 1
                    locidx += 1
                end
            elseif isa(chr, Mitochondrion)
                locidx += length(chr.loci1)
            elseif isa(chr, XChromosome) || isa(chr, YChromosome)
                for i = 1:length(chr.loci1)
                    setgene!(df, colidx, indidx, :male, d, c, locidx, chr.loci1[i])
                    locidx += 1
                    colidx += 1
                end
            end
        end
        indidx += 1
    end
    colidx, indidx
end

function getallgenes(pop, demes)
    n = nelems(pop)
    df = DataFrame(org_id = Array(Int, n),
                   sex = Array(Symbol, n),
                   deme = Array(Int8, n),
                   chr_type = Array(Symbol, n),
                   locus_id = Array(Int, n),
                   gene = Array(Gene, n))

    colidx, indidx = getgenes!(pop[1], demes[1], df, 1, 1, Female)
    getgenes!(pop[2], demes[2], df, colidx, indidx, Male)
    df
end

# function getlocidata(df)
#     loci = sort(unique(df[:locus_id]))
#     nextgid = 0
#     lociset = Array((Symbol, Int, Dict{UInt, Int}), length(loci))
#     nelems = 0
#     for (i, locus) in enumerate(loci)
#         data = df[df[:locus_id] .== locus, [:chr_type, :locus_id, :gene]]
#         c = data[1, :chr_type]
#         lociset[i] = (c, locus, Dict{UInt, Int}())
#
#         gids = sort(unique([getid(g, Mutation) for g in data[:gene]]))
#         for (j, gid) in enumerate(gids)
#             lociset[i][3][gid] = nextgid + j
#         end
#         nextgid += length(gids)
#         nelems += binomial(length(gids), 2)
#     end
#     lociset, nelems
# end

function savedistance(path, rec, df)
    loci = sort(unique(df[:locus_id]))
    nelems = 0
    gs = Array(Vector{UInt}, length(loci))
    cs = Array(ASCIIString, length(loci))
    for (i, locus) in enumerate(loci)
        subdf = df[df[:locus_id] .== locus, [:chr_type, :gene]]
        genes = sort(unique([getid(g, Mutation) for g in subdf[:gene]]))
        cs[i] = string(subdf[1, :chr_type])
        gs[i] = [int(g) for g in genes]
        nelems += binomial(length(genes), 2)
    end
    dists = DataFrame(chr_type = Array(ASCIIString, nelems),
                      locus_id = Array(Int, nelems),
                      gene_id_1 = Array(Int, nelems),
                      gene_id_2 = Array(Int, nelems),
                      distance = DataArray(Int, nelems))
    colidx = 1
    for (c, locus, g) in zip(cs, loci, gs)
        colidx = getdistance!(dists, colidx, (c, locus, g), rec)
    end
    writetable(path, dists, separator = '\t')
end

function savemigrations(path, rec, df)
    len = size(df, 1)
    gids = [getid(g, Mutation) for g in df[:gene]]
    nmigs = [countalong(rec, g, Migration) for g in df[:gene]]

    data = DataFrame(org_id = df[:org_id],
                     sex = [string(s) for s in df[:sex]],
                     deme = df[:deme],
                     chr_type = [string(c) for c in df[:chr_type]],
                     locus_id = df[:locus_id],
                     gene_id = gids,
                     number_of_migrations = nmigs)
    writetable(path, data, separator = '\t')
end

function samplemigintervals(path, rec, df)
    loci = sort(unique(df[:locus_id]))
    dfs = DataFrame[]
    for locus in loci
        gs = df[df[:locus_id] .== locus, :]
        len = size(gs, 1)
        i = rand(1:len)
        is = eventintervals(rec, gs[i, :gene], Migration)
        length(is) == 0 && push!(is, 0)
        len2 = length(is)
        sdf = DataFrame(org_id = repmat([gs[i, :org_id]], len2),
                        chr_type = repmat([gs[i, :chr_type]], len2),
                        sex = repmat([gs[i, :sex]], len2),
                        deme = repmat([gs[i, :deme]], len2),
                        locus_id = repmat([gs[i, :locus_id]], len2),
                        intervals = is)
        push!(dfs, sdf)
    end
    writetable(path, vcat(dfs...), separator = '\t')
end
