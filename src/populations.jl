using StatsBase: sample!

immutable Organisms{T<:Organism}
    size::Vector{Int}
    trait::Vector{DemeIndex}
    data::Vector{T}
end

immutable Population
    f::Organisms{Female}
    m::Organisms{Male}
end

function Population(nf::Vector{Int}, nm::Vector{Int})
    ft = vcat([fill(convert(DemeIndex, i), j) for (i, j) = enumerate(nf)]...)
    mt = vcat([fill(convert(DemeIndex, i), j) for (i, j) = enumerate(nm)]...)

    f = Organisms{Female}(nf, ft, [Female() for _ = 1:sum(nf)])
    m = Organisms{Male}(nm, mt, [Male() for _ = 1:sum(nm)])

    Population(f, m)
end

Base.getindex(p::Population, ::Type{Female}) = p.f
Base.getindex(p::Population, ::Type{Male}) = p.m

function migrate!(t, orgs, n, migs, rec)
    tar, src = listmovers(n, migs)
    b = n[1]
    for i = 1:length(tar)
        if tar[i] <= b < src[i] || src[i] <= b < tar[i]
            orgs[src[i]] = migrate!(t, orgs[src[i]], rec)
        end
    end
    orgs[tar] = orgs[src]
end

function selectmigrants!(sels, nelem, s)
    sample!(1:nelem[1], sels[1], replace = false, ordered = true)
    sample!(nelem[1]+(1:nelem[2]), sels[2], replace = false, ordered = true)
end

function selectmigrants(nelem, s)
    sels = Array(Vector{Int}, length(s))
    for i = 1:length(s)
        sels[i] = Array(Int, s[i])
    end
    selectmigrants!(sels, nelem, s)
    sels
end

function listmovers(nelem, s)
    # elems should not be used in production code, and this function will return src and tar
    sels = selectmigrants(nelem, s)

    len1 = length(sels[1])
    len2 = length(sels[2])
    if len1 == len2
        src = vcat(sels[1], sels[2])
        tar = vcat(sels[2], sels[1])
    else
        len = min(len1, len2)
        if len1 < len2
            r = nelem[1] + (1:(len2 - len))
            tar = sort(unique(vcat(sels[1], sels[2], r)))
            src = vcat(sels[2], sels[1], setdiff(r, sels[2]))
        else
            r = nelem[1] + (1-len1+len:0)
            tar = sort(unique(vcat(sels[1], sels[2], r)))
            src = vcat(setdiff(r, sels[1]), sels[2], sels[1])
        end
    end
    tar, src
end
