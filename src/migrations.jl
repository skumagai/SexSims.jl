using JSON
using DataFrames: readtable

immutable Autosome end
immutable XChromosome end
immutable YChromosome end
immutable Mitochondrion end

function m(αf, αm, βf, βm, σf, σm, τf, τm, bf, bm, chr::Type{Autosome})
    x = X(αf, αm, βf, βm, σf, σm, chr)
    first((ω(αf, αm, σf, σm, τf, x) * bf + ω(βf, βm, σf, σm, τm, x) * bm) / 2.0)
end
function m(αf, αm, βf, βm, σf, σm, τf, τm, bf, bm, chr::Type{XChromosome})
    x = X(αf, αm, βf, βm, σf, σm, chr)
    first((2.0 * ω(αf, αm, σf, σm, τf, x) * bf + ω(βf, βm, σf, σm, τm, x) * bm) / 3.0)
end
m(αf, αm, βf, βm, σf, σm, τf, τm, bf, bm, ::Type{YChromosome}) = first(ω(βm, σm, τm) * bm)
m(αf, αm, βf, βm, σf, σm, τf, τm, bf, bm, ::Type{Mitochondrion}) = first(ω(αf, σf, τf) * bf)

ω(vf, vm, σf, σm, τ, x) = τ * (1.0 - 0.5 .* [vf vm] * inv(eye(x) - x) * [1.0 - σf, 1.0 - σm])
ω(v, σ, τ) = τ * (1.0 - v) / (1.0 - σ * v)

X(af, am, bf, bm, sf, sm, ::Type{Autosome}) = diagm([sf, sm]) * diagm([0.5, 0.5]) * [af am; bf bm]
X(af, am, bf, bm, sf, sm, ::Type{XChromosome}) = diagm([sf, sm]) * diagm([0.5, 1]) * [af am; bf 0.0]

function readinput(path)
    params = JSON.parsefile(path)
    αf = float64(first(filter(x -> x["from"] == "female" && x["to"] == "female", params["learning"]))["rate"])
    αm = float64(first(filter(x -> x["from"] == "female" && x["to"] == "male", params["learning"]))["rate"])
    βf = float64(first(filter(x -> x["from"] == "male" && x["to"] == "female", params["learning"]))["rate"])
    βm = float64(first(filter(x -> x["from"] == "male" && x["to"] == "male", params["learning"]))["rate"])
    σf = float64([e[3-i] for (i, e) in enumerate(params["trait"]["female"])])
    σm = float64([e[3-i] for (i, e) in enumerate(params["trait"]["male"])])
    τf = float64(params["migration"]["female"]["cost"])
    τm = float64(params["migration"]["male"]["cost"])
    bf = float64(params["migration"]["female"]["fraction"])
    bm = float64(params["migration"]["male"]["fraction"])

    αf, αm, βf, βm, σf, σm, τf, τm, bf, bm
end

function getms(config)
    ms = Dict{ASCIIString, Vector{Float64}}(
        "Autosome" => [],
        "XChromosome" => [],
        "YChromosome" => [],
        "Mitochondrion" => []
    )
    for (αf, αm, βf, βm, σf, σm, τf, τm, bf, bm) in zip(readinput(config)...)
        push!(ms["Autosome"], m(αf, αm, βf, βm, σf, σm, τf, τm, bf, bm, Autosome))
        push!(ms["XChromosome"], m(αf, αm, βf, βm, σf, σm, τf, τm, bf, bm, XChromosome))
        push!(ms["YChromosome"], m(αf, αm, βf, βm, σf, σm, τf, τm, bf, bm, YChromosome))
        push!(ms["Mitochondrion"], m(αf, αm, βf, βm, σf, σm, τf, τm, bf, bm, Mitochondrion))
    end
    ms
end


function getnummigs(fs)
    intervals = Dict{"ASCIIString", Dict{Int, Int}}()
    for f in fs
        readtable(f)
    end
end

function main()
    length(ARGS) != 1 && error("Usage: julia migrations.jl inputfile")
    config = ARGS[1]
    isfile(config) || error("Not a file: $config")
    ms = getms(config)

    files = readdir(first(splitdir(abspath(config))))
    datadirs = filter(x -> ismatch(Regex("$config-\\d"), x) && isdir(x), files)
    datafiles = [joinpath(d, "migintervals.tsv") for d in datadirs]
    getnummigs(datafiles)
    println(ms)
end

main()
