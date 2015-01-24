# Chromosome type holds various information about chromosome such as
# the kind of chromosome as well as the number of loci on the chromosome
# and recombination rates.
# However, it doesn't hold actual genetic information.

type Chromosome
    kind::ChromosomeKind
    id::Int
    nloci::Int
    mut::Vector{Float64}
    recomb::Vector{Float64}
end

# Constructor of Chromosome type.
function Chromosome(
    k::Union(Type{Val{Autosome}},
             Type{Val{XChromosome}},
             Type{Val{YChromosome}},
             Type{Val{Mitochondrion}}),
    i::Int,
    n::Int,
    m::Vector{Float64})

    length(m) != n && error("Length of mutation rate vector is wrong.")
    m[m .< 0.0] = 0.0
    m[m .> 1.0] = 1.0

    Chromosome(k, i, n, m, fill(0.5, n - 1))
end

function Chromosome(
    k::Union(Type{Val{Autosome}},
             Type{Val{XChromosome}},
             Type{Val{YChromosome}},
             Type{Val{Mitochondrion}}),
    i::Int,
    n::Int,
    m::Vector{Float64},
    r::Vector{Float64})

    length(m) != n && error("Length of mutation rate vector is wrong.")
    length(r) != n - 1 && error("Length of recombination rate vector is wrong.")
    m[m .< 0.0] = 0.0
    r[r .> 1.0] = 1.0
    r[r .< 0.0] = 0.0
    r[r .> 0.5] = 0.5
    Chromosome(k, i, n, r)
end

# numslots return the number of slots a chromosome requires to
# store the entire genetic information.
numslots(c::Chromosome) = numslots(c.nloci, c.kind)
numslots(nloci, ::Union(Type{Val{Autosome}}, Type{Val{XChromosome}})) = 2 * nloci
numslots(nloci, ::Union(Type{Val{XChromosome}}, Type{Val{YChromosome}})) = nloci

Base.length(c::Chromosome) = length(c.nloci)

function meiosis(t, chr::Chromosome, mut, rec, ::Union(Type{Female}, Type{Male}), ::Union(Type{Female}, Type{Male}))
    l = length(chr)
    dchr = Array(Gene, l)
    for i = 1:l
        gene = rand() < 0.5 ? chr.loci1[i] : chr.loci2[i]
        dchr[i] = mutate!(t, gene, mut, rec)
    end
    tuple(dchr...)
end

meiosis(t, chr::XChromosome, mut, rec, ::Type{Male}, ::Type{Female}) = tuple([mutate!(t, gene, mut, rec) for gene in chr.loci1]...)
meiosis(t, chr::XChromosome, mut, rec, ::Type{Male}, ::Type{Male}) = chr.loci1
meiosis(t, chr::YChromosome, mut, rec, ::Type{Male}, ::Type{Male}) = tuple([mutate!(t, gene, mut, rec) for gene in chr.loci1]...)
meiosis(t, chr::YChromosome, mut, rec, ::Union(Type{Female}, Type{Male}), ::Union(Type{Female}, Type{Male})) = chr.loci1
meiosis(t, chr::Mitochondrion, mut, rec, ::Type{Female}, ::Type{Female}) = tuple([mutate!(t, gene, mut, rec) for gene in chr.loci1]...)

meiosis(t, chr::Mitochondrion, mut, rec, ::Union(Type{Female}, Type{Male}), ::Union(Type{Female}, Type{Male})) = chr.loci1
