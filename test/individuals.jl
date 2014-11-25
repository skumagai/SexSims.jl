using SexSims: GeneStateRecorder, Male, Female, Migration, Autosome, XChromosome, YChromosome, Mitochondrion, makeancestor, migrate!, getid

ret = makeancestor((Autosome{1}, XChromosome{1}, YChromosome{1}, Mitochondrion{1}))
i = 0
for (c, ex) in zip(ret, (Autosome, XChromosome, YChromosome, Mitochondrion))
    @test isa(c, ex)
    i += 1
end
@test i == 4
exit(1)

@test isa(ret.chrs[1], Autosome{1})
@test isa(ret.chrs[2], XChromosome{1})
@test isa(ret.chrs[3], YChromosome{1})
@test isa(ret.chrs[4], Mitochondrion{1})
@test length(ret) == 4

rec = GeneStateRecorder(10)
for par = (Female, Male)
    mig = migrate!(1, ret, rec, par)
    @test getid(mig.chrs[1].loci1[1], Migration) != 0x0
    @test getid(mig.chrs[1].loci2[1], Migration) != 0x0
    if par == Female
        @test getid(mig.chrs[2].loci1[1], Migration) != 0x0
        @test getid(mig.chrs[2].loci2[1], Migration) != 0x0
        @test getid(mig.chrs[3].loci1[1], Migration) == 0x0
        @test getid(mig.chrs[4].loci1[1], Migration) != 0x0
    else
        @test getid(mig.chrs[2].loci1[1], Migration) != 0x0
        @test getid(mig.chrs[2].loci2[1], Migration) == 0x0
        @test getid(mig.chrs[3].loci1[1], Migration) != 0x0
        @test getid(mig.chrs[4].loci1[1], Migration) == 0x0
    end
end
