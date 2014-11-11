using StatsBase: sample!

function main()
    nelem = [10, 10]
    elems = Int[1:20]
    s = [4, 6]
    sels = Array(Vector{Int}, length(s))
    sources = Array(Vector{Int}, length(s))
    for i = 1:length(s)
        sels[i] = Array(Int, s[i])
    end
    sources[1] = Int[nelem[1]-s[1]+(1:s[1])]
    sources[2] = Int[nelem[1]+(1:s[2])]

    sample!(sub(elems, 1:nelem[1]), sels[1], replace = false, ordered = true)
    sample!(sub(elems, nelem[1]+(1:nelem[2])), sels[2], replace = false, ordered = true)

    println(sels)
    if length(sels[1]) == length(sels[2])
        elems[[sels[1], sels[2]]] = elems[[sels[2], sels[1]]]
    else
        if length(sels[1]) < length(sels[2])
            i, j = 1, 2
            println("$i<$j")
            offset = nelem[1]
            op = +
        else
            i, j = 2, 1
            println("$i<$j")
            offset = nelem[1] + 1
            op = -
        end
        leni = length(sels[i])
        diff = op(offset, 1:(length(sels[j]) - length(sels[i])))
        src = [sels[i], sels[j][1:leni], sels[j][(leni+1):end], diff]
        tar = [sels[j][1:leni], sels[i], diff, sels[j][(leni+1):end]]
        println(src)
        println(tar)
        elems[tar] = elems[src]
    end
    println(elems)
end

main()
