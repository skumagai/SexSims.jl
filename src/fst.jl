using DataFrames: readtable, writetable, eachrow, join, DataFrame
using SexSims: readinput

function geth(data)
    d = Dict{Int, Int}()
    total = 0
    for gene in data[:gene_id]
        if gene != 0
            d[gene] = get(d, gene, 0) + 1
            total += 1
        end
    end
    total, 1 - sum([(val / total)^2 for val in values(d)])
end

function getfst(d)
    migdata = readtable(joinpath(d, "migrations.tsv"))
    loci = sort(unique(migdata[:locus_id]))
    demes = sort(unique(migdata[:deme]))
    fst1 = Array(Float64, length(loci))
    fst2 = Array(Float64, length(loci))
    fst3 = Array(Float64, length(loci))
    dname = Array(ASCIIString, length(loci))
    fill!(dname, d)
    for (i, locus) in enumerate(loci)
        ddata = [migdata[(migdata[:locus_id] .== locus) & (migdata[:deme] .== deme), :] for deme in demes]
        size1, freqs1 = geth(ddata[1])
        size2, freqs2 = geth(ddata[2])
        sizet, freqst = geth(migdata[migdata[:locus_id] .== locus, :])

        fst1[i] = (freqst - freqs1) / freqst
        fst2[i] = (freqst - freqs2) / freqst
        fst3[i] = (freqst - (size1 * freqs1 + size2 * freqs2) / (size1 + size2)) / freqst
    end
    DataFrame(dataset = dname, locus_id = loci, Fst_1 = fst1, Fst_2 = fst2, Fst_total = fst3)
end

getfsts(ds) = vcat([getfst(d) for d in ds]...)

function main()
    length(ARGS) != 1 && error("Usage: julia fst.jl inputfile")
    config = ARGS[1]
    isfile(config) || error("Note a file: $config")
    files = readdir(first(splitdir(abspath(config))))
    datadirs = filter(x -> ismatch(Regex("$config-\\d"), x) && isdir(x), files)
    fsts = getfsts(datadirs)
    writetable(joinpath(first(splitdir(abspath(config))), "Fst.tsv"), fsts, separator = '\t')
end

main()
