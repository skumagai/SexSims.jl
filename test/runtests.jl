using Base.Test

for testfile in ["types", "genes", "chromosomes", "individuals", "modelspecifics", "ios"]
    info("start test for $testfile.jl")
    include("$testfile.jl")
    info("finish test for $testfile.jl")
end
