export Gene, Change, GeneStateRecorder, id, mig, from, to, gen, migrate!, mutate!

const MUTATION = 0x0
const MIGRATION = 0x1

immutable Gene
    id::StateIndex
    mig::StateIndex
end
id(g::Gene) = g.id
mig(g::Gene) = g.mig

immutable Change
    gen::Generation
    from::StateIndex
    to::StateIndex
    kind::EventType
end
from(c::Change) = c.from
to(c::Change) = c.to
gen(c::Change) = c.gen
function kind(c::Change)
    if c.kind == MUTATION
        "mutation"
    elseif c.kind == MIGRATION
        "migration"
    else
        "undefined"
    end
end

type GeneStateRecorder
    next::StateIndex
    data::Vector{Vector{Change}}
    chunk::Int
end

function GeneStateRecorder(step)
    vec = Array(Array{Change, 1}, 1)
    vec[1] = Array(Change, step)
    GeneStateRecorder(0, vec, step)
end

function Base.getindex(st::GeneStateRecorder, i)
    j, k = divrem(i - 1, st.chunk)
    j, k = j + 1, k + 1
    st.data[j, k]
end

Base.start(st::GeneStateRecorder) = 1
function Base.next(st::GeneStateRecorder, i)
    j, k = divrem(i - 1, st.chunk)
    j, k = j + 1, k + 1
    st.data[j][k], i + 1
end
Base.done(st::GeneStateRecorder, i) = i - 1 == st.next

function prepare!(t, v, typ, st)
    st.next += 1
    i, j = divrem(st.next - 1, st.chunk)
    i, j = i + 1, j + 1
    if i > length(st.data)
        push!(st.data, Array(Change, st.chunk))
    end
    st.data[i][j] = Change(t, v, st.next, typ)
end

function mutate!(t, g::Gene, st::GeneStateRecorder)
    prepare!(t, g.id, MUTATION, st)
    Gene(st.next, g.mig)
end

function migrate!(t, g::Gene, st::GeneStateRecorder)
    prepare!(t, g.mig, MIGRATION, st)
    Gene(g.id, st.next)
end

function growifnecessary!(st::GeneStateRecorder)
    if st.next == st.cursize
        st.incr *= 2
        st.cursize += st.incr
        sizehint(st.data, st.cursize)
    end
    nothing
end
