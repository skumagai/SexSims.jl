srand(0x1)
n = 1000000
tol = 5e-2

a = makeancestor(Autosome{4})
xy = makeancestor(XYChromosome{3})
mito = makeancestor(Mitochondrion{5})
@test length(a) == 4
@test length(xy) == 3
@test length(mito) == 5

a = Autosome(((Gene(0x1, 0x0), Gene(0x2, 0x0)), (Gene(0x3, 0x0), Gene(0x4, 0x0))))
xy = XYChromosome(((Gene(0x1, 0x0), Gene(0x2, 0x0)), (Gene(0x3, 0x0), Gene(0x4, 0x0))))
mito = Mitochondrion(((Gene(0x1, 0x0), Gene(0x2, 0x0)),))

rec = GeneStateRecorder(1)
for par = (Female, Male), ch = (Female, Male)
    ms = [meiosis(1, a1, 0.0, rec, par, ch) for _ = 1:n]
    l1g1 = [m[1] for m in ms] .== Gene(0x1, 0x0)
    l1g3 = [m[1] for m in ms] .== Gene(0x3, 0x0)
    l2g2 = [m[2] for m in ms] .== Gene(0x2, 0x0)
    l2g4 = [m[2] for m in ms] .== Gene(0x4, 0x0)
    @test_approx_eq_eps countnz(l1g1) / n 0.5 tol
    @test_approx_eq_eps countnz(l1g3) / n 0.5 tol
    @test_approx_eq_eps countnz(l2g2) / n 0.5 tol
    @test_approx_eq_eps countnz(l2g4) / n 0.5 tol
    @test_approx_eq_eps countnz(l1g1 & l2g2) / n 0.25 tol
    @test_approx_eq_eps countnz(l1g3 & l2g4) / n 0.25 tol
    @test_approx_eq_eps countnz(l1g1 & l2g4) / n 0.25 tol
    @test_approx_eq_eps countnz(l1g3 & l2g2) / n 0.25 tol
end

for ch = (Female, Male)
    ms = [meiosis(1, a1, 0.0, rec, Female, ch) for _ = 1:n]
    l1g1 = [m[1] for m in ms] .== Gene(0x1, 0x0)
    l1g3 = [m[1] for m in ms] .== Gene(0x3, 0x0)
    l2g2 = [m[2] for m in ms] .== Gene(0x2, 0x0)
    l2g4 = [m[2] for m in ms] .== Gene(0x4, 0x0)
    @test_approx_eq_eps countnz(l1g1) / n 0.5 tol
    @test_approx_eq_eps countnz(l1g3) / n 0.5 tol
    @test_approx_eq_eps countnz(l2g2) / n 0.5 tol
    @test_approx_eq_eps countnz(l2g4) / n 0.5 tol
    @test_approx_eq_eps countnz(l1g1 & l2g2) / n 0.25 tol
    @test_approx_eq_eps countnz(l1g3 & l2g4) / n 0.25 tol
    @test_approx_eq_eps countnz(l1g1 & l2g4) / n 0.25 tol
    @test_approx_eq_eps countnz(l1g3 & l2g2) / n 0.25 tol
end

ms = [meiosis(1, xy1, 0.0, rec, Male, Female) for _ = 1:n]
l1g1 = [m[1] for m in ms] .== Gene(0x1, 0x0)
l2g2 = [m[2] for m in ms] .== Gene(0x2, 0x0)
@test countnz(l1g1) / n == 1.0
@test countnz(l2g2) / n == 1.0
@test countnz(l1g1 & l2g2) / n == 1.0

ms = [meiosis(1, xy1, 0.0, rec, Male, Male) for _ = 1:n]
l1g3 = [m[1] for m in ms] .== Gene(0x3, 0x0)
l2g4 = [m[2] for m in ms] .== Gene(0x4, 0x0)
@test countnz(l1g3) / n == 1.0
@test countnz(l2g4) / n == 1.0
@test countnz(l1g3 & l2g4) / n == 1.0

ms = [meiosis(1, mito, 0.5, rec, Female, Female) for _ = 1:n]
l1g1 = [m[1] for m in ms] .== Gene(0x1, 0x0)
l2g2 = [m[2] for m in ms] .== Gene(0x2, 0x0)
@test_approx_eq_eps countnz(l1g1) / n 0.5 tol
@test_approx_eq_eps countnz(l2g2) / n 0.5 tol
@test_approx_eq_eps countnz(l1g1 & l2g2) / n 0.25 tol

ms = [meiosis(1, mito, 0.5, rec, Female, Male) for _ = 1:n]
l1g1 = [m[1] for m in ms] .== Gene(0x1, 0x0)
l2g2 = [m[2] for m in ms] .== Gene(0x2, 0x0)
@test countnz(l1g1) / n == 1.0
@test countnz(l2g2) / n == 1.0
@test countnz(l1g1 & l2g2) / n == 1.0

for ch = (Female, Male)
    ms = [meiosis(1, mito, 0.0, rec, Male, ch) for _ = 1:n]
    l1g1 = [m[1] for m in ms] .== Gene(0x1, 0x0)
    l2g2 = [m[2] for m in ms] .== Gene(0x2, 0x0)
    @test countnz(l1g1) / n == 1.0
    @test countnz(l2g2) / n == 1.0
    @test countnz(l1g1 & l2g2) / n == 1.0
end

@test fuse(a.loci[1], a.loci[2], typeof(a)) == a
@test fuse(xy.loci[1], xy.loci[2], typeof(xy)) == xy
@test fuse(mito.loci[1], mito.loci[1], typeof(mito)) == mito

for chr = (Autosome, XYChromosome), par = (Female, Male)
    ret = migrate!(1, makeancestor(chr{1}), rec, par)
    @test SexSims.getid(ret.loci[1][1], SexSims.Migration) != 0x0
    @test SexSims.getid(ret.loci[2][1], SexSims.Migration) != 0x0
end

ret = migrate!(1, makeancestor(Mitochondrion{1}), rec, Female)
@test SexSims.getid(ret.loci[1][1], SexSims.Migration) != 0x0
ret = migrate!(1, makeancestor(Mitochondrion{1}), rec, Male)
@test SexSims.getid(ret.loci[1][1], SexSims.Migration) == 0x0
