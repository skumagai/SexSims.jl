using SexSims: GeneStateRecorder, Male, Female, Migration, Autosome, XYChromosome, Mitochondrion, makeancestor, migrate!, getid

ret = makeancestor((Autosome{1}, XYChromosome{1}, Mitochondrion{1}))
@test isa(ret.chrs[1], Autosome{1})
@test isa(ret.chrs[2], XYChromosome{1})
@test isa(ret.chrs[3], Mitochondrion{1})
@test length(ret) == 3

rec = GeneStateRecorder(10)
for par = (Female, Male)
    mig = migrate!(1, ret, rec, par)
    for i = 1:2
        @test getid(mig.chrs[i].loci1[1], Migration) != 0x0
        @test getid(mig.chrs[i].loci2[1], Migration) != 0x0
    end
    if par == Female
        @test getid(mig.chrs[3].loci1[1], Migration) != 0x0
    else
        @test getid(mig.chrs[3].loci1[1], Migration) == 0x0
    end
end
