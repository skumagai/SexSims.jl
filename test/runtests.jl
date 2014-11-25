using Base.Test

# for testfile in ["types", "genes", "chromosomes", "individuals", "modelspecifics", "ios"]
for testfile in ["genes"]
    info("start test for $testfile.jl")
    include("$testfile.jl")
    info("finish test for $testfile.jl")
end
