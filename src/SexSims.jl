module SexSims

export Female, Male, Population, ParentIds, GeneStateRecorder, migrate!, fitness!, learn!, reproduce!, pickparents!

include("types.jl")
include("genes.jl")
include("individuals.jl")
include("populations.jl")

end
