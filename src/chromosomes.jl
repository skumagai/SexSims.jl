export Autosome, XYChromosome, Mitochondrion
export makeancestor, meiosis, fuse

typealias Loci{N} NTuple{N, Gene}

abstract Chromosome

immutable Autosome{N} <: Chromosome
    loci::NTuple{2, Loci{N}}
end

immutable XYChromosome{N} <: Chromosome
    loci::NTuple{2, Loci{N}}
end

immutable Mitochondrion{N} <: Chromosome
    loci::NTuple{1, Loci{N}}
end

function makeancestor{T<:Chromosome}(chr::Type{T})
    chr1 = [Gene(0x0, 0x0) for _ = 1:chr.parameters[1]]
    chr2 = copy(chr1)
    chr((tuple(chr1...), tuple(chr2...)))
end

function makeancestor{T<:Mitochondrion}(chr::Type{T})
    chr((tuple([Gene(0x0, 0x0) for i = 1:chr.parameters[1]]...),))
end

Base.length(c::Chromosome) = length(c.loci[1])

function meiosis{S1<:Sex, S2<:Sex}(t, chr::Autosome, mut, rec, ::Type{S1}, ::Type{S2})
    dchr = Array(Gene, length(chr))
    for (i, (gene1, gene2)) = enumerate(zip(chr.loci[1], chr.loci[2]))
        gene = rand() < 0.5 ? gene1 : gene2
        dchr[i] = mutate!(t, gene, mut, rec)
    end
    tuple(dchr...)
end

function meiosis(t, chr::XYChromosome, mut, rec, ::Type{Male}, ::Type{Female})
    tuple([mutate!(t, gene, mut, rec) for gene in chr.loci[1]]...)
end

function meiosis(t, chr::XYChromosome, mut, rec, ::Type{Male}, ::Type{Male})
    tuple([mutate!(t, gene, mut, rec) for gene in chr.loci[2]]...)
end

function meiosis(t, chr::Mitochondrion, mut, rec, ::Type{Female}, ::Type{Female})
    tuple([mutate!(t, gene, mut, rec) for gene in chr.loci[1]]...)
end

function meiosis{S1<:Sex, S2<:Sex}(t, chr::Mitochondrion, mut, rec, ::Type{S1}, ::Type{S2})
    chr.loci[1]
end

fuse{T<:Chromosome}(chr1, chr2, ::Type{T}) = T((chr1, chr2))
fuse{T<:Mitochondrion}(chr1, chr2, ::Type{T}) = T((chr1,))

function migrate!{S<:Sex}(t, chr::Chromosome, rec, sex::Type{S})
    typeof(chr)((migrate!(t, chr.loci[1], rec, sex), migrate!(t, chr.loci[2], rec, sex)))
end

function migrate!(t, chr::Mitochondrion, rec, ::Type{Female})
    Mitochondrion((migrate!(t, chr.loci[1], rec, Female),))
end

migrate!(t, chr::Mitochondrion, rec, ::Type{Male}) = chr

function migrate!{S<:Sex}(t, loci::(Gene...), rec, ::Type{S})
    tuple([migrate!(t, locus, rec) for locus in loci]...)
end
