using JSON
using PyPlot
using DataFrames: readtable
using Distributions: Geometric, pdf
using LaTeXStrings
using PyPlot
using PyCall

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
    Nf = params["population size"]["female"]
    Nm = params["population size"]["male"]
    Dict{UTF8String, Any}(
        "αf" => αf,
        "αm" => αm,
        "βf" => βf,
        "βm" => βm,
        "σf" => σf,
        "σm" => σm,
        "τf" => τf,
        "τm" => τm,
        "bf" => bf,
        "bm" => bm,
        "Nf" => Nf,
        "Nm" => Nm
    )
end

function getms(p)
    ms = Dict{ASCIIString, Vector{Float64}}(
        "autosome" => [],
        "xchromosome" => [],
        "ychromosome" => [],
        "mitochondrion" => []
    )
    for i = 1:2
        push!(ms["autosome"], m(p["αf"][i], p["αm"][i], p["βf"][i], p["βm"][i], p["σf"][i], p["σm"][i], p["τf"][i], p["τm"][i], p["bf"][i], p["bm"][i], Autosome))
        push!(ms["xchromosome"], m(p["αf"][i], p["αm"][i], p["βf"][i], p["βm"][i], p["σf"][i], p["σm"][i], p["τf"][i], p["τm"][i], p["bf"][i], p["bm"][i], XChromosome))
        push!(ms["ychromosome"], m(p["αf"][i], p["αm"][i], p["βf"][i], p["βm"][i], p["σf"][i], p["σm"][i], p["τf"][i], p["τm"][i], p["bf"][i], p["bm"][i], YChromosome))
        push!(ms["mitochondrion"], m(p["αf"][i], p["αm"][i], p["βf"][i], p["βm"][i], p["σf"][i], p["σm"][i], p["τf"][i], p["τm"][i], p["bf"][i], p["bm"][i], Mitochondrion))
    end
    ms
end

function getnummigs(fs)
    intervals = Dict{Int, (ASCIIString, NTuple{2, Array{Int, 1}})}()
    df = readtable(fs[1])
    loci = unique(df[:locus_id])
    for locus in loci
        deme = first(df[df[:locus_id] .== locus, :deme])
        chr = first(df[df[:locus_id] .== locus, :chr_type])
        vals = df[df[:locus_id] .== locus, :intervals]
        intervals[locus] = (chr, (vals[deme:2:end], vals[(3-deme):2:end]))
    end
    for f in fs[2:end]
        df = readtable(f)
        loci = unique(df[:locus_id])
        for locus in loci
            deme = first(df[df[:locus_id] .== locus, :deme])
            chr = first(df[df[:locus_id] .== locus, :chr_type])
            vals = df[df[:locus_id] .== locus, :intervals]
            append!(intervals[locus][2][1], vals[deme:2:end])
            append!(intervals[locus][2][2], vals[(3-deme):2:end])
        end
    end
    intervals
end

function plotresults(p, model, data, fname)
    nloci = length(data)
    PyPlot.figure(figsize = (9, 14))
    # ymax = 0.0
    # ymax2 = 0.0
    for (i, d) in enumerate(keys(data))
        chr = data[d][1]
        for (j, vals) in enumerate(data[d][2])
            xmax = maximum(vals)
            PyPlot.subplot(nloci + 1, 2, (i - 1) * 2 + j)
            PyPlot.plt.hist(vals, bins = 20, normed = true)
            # ymax2 = PyPlot.plt.ylim()[2]
            # ymax = ymax > ymax2 ? ymax : ymax2

            g = Geometric(model[data[d][1]][j])
            xs = [0:xmax]
            PyPlot.plt.plot(xs, [pdf(g, x) for x in xs])
            PyPlot.plt.legend(["analytical", "simulation"])
            PyPlot.plt.title("$chr, deme = $j")
            # ymax2 = PyPlot.plt.ylim()[2]
            # ymax = ymax > ymax2 ? ymax : ymax2
        end
    end
    # for i = 1:nloci, j = 1:2
    #     PyPlot.plt.subplot(nloci + 1, 2, (i - 1) * 2 + j)
    #     PyPlot.plt.ylim(0.0, ymax)
    # end
    for i = 1:2
        ax = PyPlot.subplot(nloci + 1, 2, 2nloci + i)
        PyPlot.text(0.1, 0.95, latexstring("\$N_{f,$i} = $(p["Nf"][i])\$"))
        PyPlot.text(0.6, 0.95, latexstring("\$N_{m,$i} = $(p["Nm"][i])\$"))
        PyPlot.text(0.1, 0.85, latexstring("\$\\alpha_{f,$i} = $(p["αf"][i])\$"))
        PyPlot.text(0.6, 0.85, latexstring("\$\\alpha_{m,$i} = $(p["αm"][i])\$"))
        PyPlot.text(0.1, 0.75, latexstring("\$\\beta_{f,$i} = $(p["βf"][i])\$"))
        PyPlot.text(0.6, 0.75, latexstring("\$\\beta_{m,$i} = $(p["βm"][i])\$"))
        PyPlot.text(0.1, 0.65, latexstring("\$\\sigma_{f,$i} = $(p["σf"][i])\$"))
        PyPlot.text(0.6, 0.65, latexstring("\$\\sigma_{m,$i} = $(p["σm"][i])\$"))
        PyPlot.text(0.1, 0.55, latexstring("\$\\tau_{f,$i} = $(p["τf"][i])\$"))
        PyPlot.text(0.6, 0.55, latexstring("\$\\tau_{m,$i} = $(p["τm"][i])\$"))
        PyPlot.text(0.1, 0.45, latexstring("\$b_{f,$i} = $(p["bf"][i])\$"))
        PyPlot.text(0.6, 0.45, latexstring("\$b_{m,$i} = $(p["bm"][i])\$"))

        PyPlot.text(0.25, 0.35, latexstring("\$m^A_$i = $(@sprintf("%.3e", model["autosome"][i]))\$"))
        PyPlot.text(0.25, 0.25, latexstring("\$m^X_$i = $(@sprintf("%.3e", model["xchromosome"][i]))\$"))
        PyPlot.text(0.25, 0.15, latexstring("\$m^Y_$i = $(@sprintf("%.3e", model["ychromosome"][i]))\$"))
        PyPlot.text(0.25, 0.05, latexstring("\$m^M_$i = $(@sprintf("%.3e", model["mitochondrion"][i]))\$"))

        ax[:set_frame_on](false)
        ax[:axes][:get_xaxis]()[:set_visible](false)
        ax[:axes][:get_yaxis]()[:set_visible](false)
    end
    PyPlot.plt.savefig(fname)
end

function main()
    length(ARGS) != 1 && error("Usage: julia migrations.jl inputfile")
    config = ARGS[1]
    isfile(config) || error("Not a file: $config")
    params = readinput(config)
    ms = getms(params)

    files = readdir(first(splitdir(abspath(config))))
    datadirs = filter(x -> ismatch(Regex("$config-\\d"), x) && isdir(x), files)
    datafiles = [joinpath(d, "migintervals.tsv") for d in datadirs]
    migs = getnummigs(datafiles)
    plotresults(params, ms, migs, splitext(config)[1]".pdf")
end

main()
