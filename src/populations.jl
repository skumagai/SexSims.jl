# Population holds a vector (1D array) of genes.
# A individual is not explicitly modeled, but it
# is represented as a contigious block of genes within
# the vector.
# A gene on autosome and x-chromosome uses two consective
# slots to represent genotype, whereas other types of
# chromosomes (e.g. y-chromosome and mitochondrion) only
# uses one slot.
#
# In addition to genetic data, Population type holds two
# additional information.  A vector of objects of Chromosome
# type, and the number of individuals.
type Population
    data::Vector{Gene}
    chroms::Vector{Chromosome}
    indsize::Int
    n::Int
end

# Constructor of Population type.  All genes are initialized as
# having ancestoral trait (0) as well as no migration so far.
function Population(n::Integer, chroms::Chromosome...)
    ngenes = 0
    for chrom in chroms
        ngenes += numslots(chrom)
    end
    d = fill(Gene(), n * ngenes)
    cs = [c for c in chroms]
    Population(d, cs, ngenes, n)
end

Base.length(p::Population) = p.n

# Return the first gene of an individual with respect to the entire population.
getindividualoffset(p::Population, i::Int) = (i - 1) * p.indsize + 1

function migrate!(
    t::Int,
    rec::GeneStateRecorder,
    p::Population,
    i::Int,
    sex::Sex)

    i = getindividualoffset(p, i)
    for c in p.chroms, _ = 1:length(c)
        i = migrate!(t, rec, Val{c.kind}, p, i, Val{sex})
    end
end

function migrate!(
    t::Int,
    rec::GeneStateRecorder,
    ::Type{Val{Autosome}},
    p::Population,
    i::Int,
    ::Union(Type{Val{Female}},
            Type{Val{Male}}))

    for j = i:(i+1)
        p[j] = migrate!(t, rec, p[j])
    end
    j + 1
end

function migrate!(
    t::Int,
    rec::GeneStateRecorder,
    ::Type{Val{XChromosome}},
    p::Population,
    i::Int,
    ::Type{Val{Female}})

    for j = i:(i+1)
        p[j] = migrate!(t, rec, p[j])
    end
    j + 1
end

function migrate!(
    t::Int,
    rec::GeneStateRecorder,
    ::Type{Val{XChromosome}},
    p::Population,
    i::Int,
    ::Type{Val{Male}})

    p[i] = migrate!(t, rec, p[i])
    i + 2
end

function migrate!(
    t::Int,
    rec::GeneStateRecorder,
    ::Type{Val{YChromosome}},
    p::Population,
    i::Int,
    ::Type{Val{Female}})

    i + 1
end

function migrate!(
    t::Int,
    rec::GeneStateRecorder,
    ::Type{Val{YChromosome}},
    p::Population,
    i::Int,
    ::Type{Val{Male}})

    p[i] = migrate!(t, rec, p[i])
    i + 1
end

function migrate!(
    t::Int,
    rec::GeneStateRecorder,
    ::Type{Val{Mitochondrion}},
    p::Population,
    i::Int,
    ::Type{Val{Female}})

    p[i] = migrate!(t, rec, p[i])
    i + 1
end

function migrate!(
    t::Int,
    rec::GeneStateRecorder,
    ::Type{Val{YChromosome}},
    p::Population,
    i::Int,
    ::Type{Val{Male}})

    i + 1
end

function meiosis!(
    t::Int,
    rec::GeneStateRecorder,
    pp::Populatoin,
    ip::Int,
    psex::Sex,
    pc::Population,
    ic::Int,
    csex::Sex)

    # Get offset of parent and offspring individuals
    ip = getindividualoffset(pp, ip)
    ic = getindividualoffset(pc, ic)
    for c in p.chroms
        ip, ic = meiosis!(t, rec, Val{c.kind}, pp, ip, Val{psex}, pc, ic, Val{csex})
    end
end

function meiosis!(
    t::Int,
    rec::GeneStateRecorder,
    ::Type{Val{Autosome}},
    pp::Population,
    ip::Int,
    ::Union(::Type{Val{Female}}, ::Type{Val{Male}}),
    pc::Population,
    ic::Int,
    ::Type{Val{Female}})


