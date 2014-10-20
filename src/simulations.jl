using StatsBase: sample

function parentids!(p::Vector{(DemeIndex, Int)}, v1::WeightVec, v2::WeightVec, frac::Float64)
    for i = 1:length(v1)
        if rand() < frac
            p[i] = (1, sample(v1))
        else
            p[i] = (2, sample(v2))
        end
    end
end

function parentids!(p::(Vector{(DemeIndex, Int)}, Vector{(DemeIndex, Int)}), f::(Vector{Float64}, Vector{Float64}), b::(Float64, Float64))
    parentids!(p[1],
end

function learn!(deme::DemeIndex, ctrait::Vector{DemeIndex}, momid::Vector{Int}, popid::Vector{Int}, momtrait::Vector{DemeIndex}, poptrait::Vector{DemeIndex}, momrate::Float64, poprate::Float64)
    for i = 1:length(momid)
        f = momid[i]
        m = popid[i]
        ctrait[i] = learn(deme, momtrait[f], poptrait[m], momrate, poprate)
    end
    nothing
end

function simulation()
    srand(1234)
    nf = (100, 200)
    nm = (200, 100)
    bf = (0.2, 0.3)
    bm = (0.3, 0.2)

    pop1 = Population(nf, nm)
    pop2 = Population(nf, nm)
    rec = GeneStateRecorder(4 * (sum(nf) + sum(nm)))

    # Index of parents
    momtof = (Array((DemeIndex, Int), nf[1]), Array((DemeIndex, Int), nf[2]))
    poptof = (Array((DemeIndex, Int), nf[1]), Array((DemeIndex, Int), nf[2]))
    momtom = (Array((DemeIndex, Int), nm[1]), Array((DemeIndex, Int), nm[2]))
    poptof = (Array((DemeIndex, Int), nm[1]), Array((DemeIndex, Int), nm[2]))
    for t = 1:tmax
        # mating
        parentids!(momtof, ffit, bf)
        parentids!(poptof, mfit, bf)
        parentids!(momtom, ffit, bm)
        parentids!(poptom, mfit, bm)

        parentids!(momtom1, femtrait1, femtrait2, bf1)
        parentids!(momtom2, femtrait1, femmtrait2, 1.0 - bf2)
        parentids!(poptom1, maletrait1, maletrait2, bm1)
        parentids!(poptom2, maletrait1, maletrait2, 1.0 - bm2)
        # learning

        # migrating
    end
end
