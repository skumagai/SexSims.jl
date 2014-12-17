immutable EventType
    event::UInt8
    str::ASCIIString
end
const Mutation = EventType(0x1, "mutation")
const Migration = EventType(0x2, "migration")

immutable Gene
    id::UInt
    mig::UInt
end

Gene() = Gene(uint(0), uint(0))

function getid(g::Gene, event)
    if event == Mutation
        g.id
    elseif event == Migration
        g.mig
    else
        uint(0)
    end
end

immutable Change
    gen::Int
    from::UInt
    kind::EventType
end
parent(c::Change) = c.from
generation(c::Change) = c.gen
kind(c::Change) = c.kind.str

type GeneStateRecorder
    next::UInt
    data::Vector{Vector{Change}}
    chunk::UInt
end

function GeneStateRecorder(step)
    vec = Array(Vector{Change}, 1)
    vec[1] = Array(Change, step)
    GeneStateRecorder(0, vec, step)
end

function Base.getindex(st::GeneStateRecorder, i)
    j, k = divrem(i - 1, st.chunk)
    j, k = j + 1, k + 1
    st.data[j][k]
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
    st.data[i][j] = Change(t, v, typ)
end

mutate!(t, g, mut, st) = rand() < mut ? mutate!(t, g, st) : g

function mutate!(t, g::Gene, st::GeneStateRecorder)
    prepare!(t, g.id, Mutation, st)
    Gene(st.next, g.mig)
end

function migrate!(t, g::Gene, st::GeneStateRecorder)
    prepare!(t, g.mig, Migration, st)
    Gene(g.id, st.next)
end

countalong(st, gene, event) = countalong(st, getid(gene, event))

function countalong(st, idx)
    n = 0
    while idx != ANCESTRAL
        idx = st[idx].from
        n += 1
    end
    n
end

 eventintervals(st, gene, event) = eventintervals(st, getid(gene, event))

function eventintervals(st, idx)
    n = countalong(st, idx)
    if n > 2
        data = Array(Int, n - 1)
        next = st[idx].from
        i = 1
        while next != ANCESTRAL
            data[i] = generation(st[idx]) - generation(st[next])
            i += 1
            next, idx = st[next].from, next
        end
        data
    else
        Int[]
    end
end

listevents(st, gene, event) = listevents(st, getid(gene, event))

function listevents(st, idx)
    n = countalong(st, idx)
    data = Array(UInt, n)
    i = 1
    while idx != ANCESTRAL
        data[i] = idx
        idx = st[idx].from
        i += 1
    end
    data
end

distance(st, gene1, gene2, event) = distance(st, getid(gene1, event), getid(gene2, event))

function distance(st, idx1, idx2)
    data1 = listevents(st, idx1)
    data2 = listevents(st, idx2)
    isc = intersect(data1, data2)
    if isempty(isc)
        Nullable{Int}()
    else
        i = isc[1]
        Nullable(findfirst(data1, i) + findfirst(data2, i) - 1)
    end
end
