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
    data::Vector{Change}
    incr::Int
    cursize::Int
end

function GeneStateRecorder(step)
    vec = Array(Change, 0)
    sizehint(vec, step)
    GeneStateRecorder(0, vec, step, step)
end

function mutate!(t, orig::Gene, st::GeneStateRecorder)
    growifnecessary!(st)
    st.next += 1
    push!(st.data, Change(t, orig.id, st.next, MUTATION))
    Gene(st.next, orig.mig)
end

function migrate!(t, g::Gene, st::GeneStateRecorder)
    growifnecessary!(st)
    st.next += 1
    push!(st.data, Change(t, g.mig, st.next, MIGRATION))
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
