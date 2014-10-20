const MUTATION = 0x0
const MIGRATION = 0x1

immutable Gene
    id::StateIndex
    mig::StateIndex
end
id(g::Gene) = g.id
mig(g::Gene) = g.mig

immutable Change
    from::StateIndex
    to::StateIndex
    gen::Generation
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

function GeneStateRecorder(step::Int)
    vec = Array(Change, 0)
    sizehint(vec, step)
    GeneStateRecorder(0, vec, step, step)
end

function mutate!(t::Generation, orig::Gene, st::GeneStateRecorder)
    growifnecessary!(st)
    st.next += 1
    push!(st.data, Change(orig.id, st.next, t, MUTATION))
    Gene(st.next, orig.mig)
end

function migrate!(t::Generation, g::Gene, st::GeneStateRecorder)
    growifnecessary!(st)
    st.next += 1
    push!(st.data, Change(g.mig, st.next, MIGRATION))
    Gene(g.id, st.next)
end

function growifnecessary!(st::GeneStateRecorder)
    if length(st.data) == st.cursize
        st.incr *= 2
        st.cursize += st.incr
        sizehint(st.data, st.cursize)
    end
    nothing
end
