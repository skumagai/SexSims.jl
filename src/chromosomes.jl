typealias Loci{N} NTuple{N, Gene}

abstract Chromosome

immutable Autosome{N} <: Chromosome
    loci1::Loci{N}
    loci2::Loci{N}
end

immutable XChromosome{N} <: Chromosome
    loci1::Loci{N}
    loci2::Loci{N}
end

immutable YChromosome{N} <: Chromosome
    loci1::Loci{N}
end
YChromosome{N}(loc1::Loci{N}, loc2::Loci{N}) = YChromosome{N}(loc2)

immutable Mitochondrion{N} <: Chromosome
    loci1::Loci{N}
end
Mitochondrion{N}(loc1::Loci{N}, loc2::Loci{N}) = Mitochondrion{N}(loc1)

function makeancestor{T<:Chromosome}(chr::Type{T})
    chr1 = [Gene() for _ = 1:chr.parameters[1]]
    chr2 = copy(chr1)
    chr(tuple(chr1...), tuple(chr2...))
end

makeancestor{T<:Union(YChromosome, Mitochondrion)}(chr::Type{T}) = chr(tuple([Gene() for _ = 1:chr.parameters[1]]...))

Base.length(c::Chromosome) = length(c.loci1)

function meiosis(t, chr::Chromosome, mut, rec, ::Union(Type{Female}, Type{Male}), ::Union(Type{Female}, Type{Male}))
    dchr = Array(Gene, length(chr))
    for (i, (gene1, gene2)) = enumerate(zip(chr.loci1, chr.loci2))
        gene = rand() < 0.5 ? gene1 : gene2
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

fuse{T<:Chromosome}(chr1, chr2, ::Type{T}) = T(chr1, chr2)
fuse{T<:YChromosome}(chr1, chr2, ::Type{T}) = T(chr2)
fuse{T<:Mitochondrion}(chr1, chr2, ::Type{T}) = T(chr1)

migrate!(t, chr::Chromosome, rec, ::Union(Type{Female}, Type{Male})) = typeof(chr)(migrate!(t, chr.loci1, rec), migrate!(t, chr.loci2, rec))
migrate!(t, chr::XChromosome, rec, ::Type{Male}) = typeof(chr)(migrate!(t, chr.loci1, rec), chr.loci2)
migrate!(t, chr::YChromosome, rec, ::Type{Female}) = chr
migrate!(t, chr::YChromosome, rec, ::Type{Male}) = typeof(chr)(migrate!(t, chr.loci1, rec))
migrate!(t, chr::Mitochondrion, rec, ::Type{Female}) = typeof(chr)(migrate!(t, chr.loci1, rec))
migrate!(t, chr::Mitochondrion, rec, ::Type{Male}) = chr
migrate!(t, loci::(Gene...), rec) = tuple([migrate!(t, locus, rec) for locus in loci]...)

getgenotype(chr::Chromosome) = [(loc1, loc2) for (loc1, loc2) in zip(chr.loci1, chr.loci2)]
getgenotype(chr::Union(YChromosome, Mitochondrion)) = [(loc1,) for loc1 in chr.loci1]
