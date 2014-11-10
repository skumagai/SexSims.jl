using StatsBase: WeightVec

immutable Organisms{T<:Organism}
    size::(Int, Int)
    trait::Vector{DemeIndex}
    data::Vector{T}
end

immutable Population
    f::Organisms{Female}
    m::Organisms{Male}
end

function Population(nf::(Int, Int), nm::(Int, Int))
    ft = vcat([fill(convert(DemeIndex, i), j) for (i, j) = enumerate(nf)]...)
    mt = vcat([fill(convert(DemeIndex, i), j) for (i, j) = enumerate(nm)]...)

    f = Organisms{Female}(nf, ft, [Female() for _ = 1:sum(nf)])
    m = Organisms{Male}(nm, mt, [Male() for _ = 1:sum(nm)])

    Population(f, m)
end

Base.getindex(p::Population, ::Type{Female}) = p.f
Base.getindex(p::Population, ::Type{Male}) = p.m
