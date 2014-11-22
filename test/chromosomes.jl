using SexSims: Gene, GeneStateRecorder, Female, Male, Autosome, XChromosome, YChromosome, Mitochondrion, Migration, makeancestor, meiosis, fuse, migrate!, getgenotype, getid
srand(0x1)
n = 100000
tol = 5e-2

a = makeancestor(Autosome{4})
x = makeancestor(XChromosome{3})
y = makeancestor(YChromosome{6})
mito = makeancestor(Mitochondrion{5})
@test length(a) == 4
@test length(x) == 3
@test length(y) == 6
@test length(mito) == 5

a = Autosome((Gene(0x1, 0x0), Gene(0x2, 0x0)), (Gene(0x3, 0x0), Gene(0x4, 0x0)))
x = XChromosome((Gene(0x1, 0x0), Gene(0x2, 0x0)), (Gene(0x3, 0x0), Gene(0x4, 0x0)))
y = YChromosome((Gene(0x1, 0x0), Gene(0x2, 0x0)))
y2 = YChromosome((Gene(0x3, 0x0), Gene(0x4, 0x0)), (Gene(0x5, 0x0), Gene(0x6, 0x0)))
mito = Mitochondrion((Gene(0x1, 0x0), Gene(0x2, 0x0)),)
mito2 = Mitochondrion((Gene(0x3, 0x0), Gene(0x4, 0x0)), (Gene(0x5, 0x0), Gene(0x6, 0x0)))

rec = GeneStateRecorder(5n)
# for autosomes
for par = (Female, Male), ch = (Female, Male)
    ms = [meiosis(1, a, 0.0, rec, par, ch) for _ = 1:n]
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

    ms = [meiosis(1, a, 0.5, rec, par, ch) for _ = 1:n]
    l1g1 = [m[1] for m in ms] .== Gene(0x1, 0x0)
    l1g3 = [m[1] for m in ms] .== Gene(0x3, 0x0)
    l2g2 = [m[2] for m in ms] .== Gene(0x2, 0x0)
    l2g4 = [m[2] for m in ms] .== Gene(0x4, 0x0)
    @test_approx_eq_eps countnz(l1g1) / n 0.25 tol
    @test_approx_eq_eps countnz(l1g3) / n 0.25 tol
    @test_approx_eq_eps countnz(l2g2) / n 0.25 tol
    @test_approx_eq_eps countnz(l2g4) / n 0.25 tol
    @test_approx_eq_eps countnz(l1g1 & l2g2) / n 0.0625 tol
    @test_approx_eq_eps countnz(l1g3 & l2g4) / n 0.0625 tol
    @test_approx_eq_eps countnz(l1g1 & l2g4) / n 0.0625 tol
    @test_approx_eq_eps countnz(l1g3 & l2g2) / n 0.0625 tol
end

# for X-chromosome
for ch = (Female, Male)
    ms = [meiosis(1, x, 0.0, rec, Female, ch) for _ = 1:n]
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

    ms = [meiosis(1, x, 0.5, rec, Female, ch) for _ = 1:n]
    l1g1 = [m[1] for m in ms] .== Gene(0x1, 0x0)
    l1g3 = [m[1] for m in ms] .== Gene(0x3, 0x0)
    l2g2 = [m[2] for m in ms] .== Gene(0x2, 0x0)
    l2g4 = [m[2] for m in ms] .== Gene(0x4, 0x0)
    @test_approx_eq_eps countnz(l1g1) / n 0.25 tol
    @test_approx_eq_eps countnz(l1g3) / n 0.25 tol
    @test_approx_eq_eps countnz(l2g2) / n 0.25 tol
    @test_approx_eq_eps countnz(l2g4) / n 0.25 tol
    @test_approx_eq_eps countnz(l1g1 & l2g2) / n 0.0625 tol
    @test_approx_eq_eps countnz(l1g3 & l2g4) / n 0.0625 tol
    @test_approx_eq_eps countnz(l1g1 & l2g4) / n 0.0625 tol
    @test_approx_eq_eps countnz(l1g3 & l2g2) / n 0.0625 tol

    ms = [meiosis(1, x, 0.0, rec, Male, ch) for _ = 1:n]
    l1g1 = [m[1] for m in ms] .== Gene(0x1, 0x0)
    l2g2 = [m[2] for m in ms] .== Gene(0x2, 0x0)
    @test countnz(l1g1) / n == 1.0
    @test countnz(l2g2) / n == 1.0
    @test countnz(l1g1 & l2g2) / n == 1.0
end

ms = [meiosis(1, x, 0.5, rec, Male, Male) for _ = 1:n]
l1g1 = [m[1] for m in ms] .== Gene(0x1, 0x0)
l2g2 = [m[2] for m in ms] .== Gene(0x2, 0x0)
@test countnz(l1g1) / n == 1.0
@test countnz(l2g2) / n == 1.0
@test countnz(l1g1 & l2g2) / n == 1.0

# for Y-chromosome
for ch = (Female, Male)
    for mut = [0.0, 0.5]
        ms = [meiosis(1, y, mut, rec, Female, ch) for _ = 1:n]
        l1g1 = [m[1] for m in ms] .== Gene(0x1, 0x0)
        l2g2 = [m[2] for m in ms] .== Gene(0x2, 0x0)
        @test countnz(l1g1) / n == 1.0
        @test countnz(l2g2) / n == 1.0
        @test countnz(l1g1 & l2g2) / n == 1.0
    end

    ms = [meiosis(1, y, 0.0, rec, Male, ch) for _ = 1:n]
    l1g1 = [m[1] for m in ms] .== Gene(0x1, 0x0)
    l2g2 = [m[2] for m in ms] .== Gene(0x2, 0x0)
    @test countnz(l1g1) / n == 1.0
    @test countnz(l2g2) / n == 1.0
    @test countnz(l1g1 & l2g2) / n == 1.0
end

ms = [meiosis(1, y, 0.5, rec, Male, Female) for _ = 1:n]
l1g1 = [m[1] for m in ms] .== Gene(0x1, 0x0)
l2g2 = [m[2] for m in ms] .== Gene(0x2, 0x0)
@test countnz(l1g1) / n == 1.0
@test countnz(l2g2) / n == 1.0
@test countnz(l1g1 & l2g2) / n == 1.0

ms = [meiosis(1, y, 0.5, rec, Male, Male) for _ = 1:n]
l1g1 = [m[1] for m in ms] .== Gene(0x1, 0x0)
l2g2 = [m[2] for m in ms] .== Gene(0x2, 0x0)
@test_approx_eq_eps countnz(l1g1) / n 0.5 tol
@test_approx_eq_eps countnz(l2g2) / n 0.5 tol
@test_approx_eq_eps countnz(l1g1 & l2g2) / n 0.25 tol

# for mitochondrioal chromosome
for ch = (Female, Male)
    for mut = [0.0, 0.5]
        ms = [meiosis(1, mito, mut, rec, Male, ch) for _ = 1:n]
        l1g1 = [m[1] for m in ms] .== Gene(0x1, 0x0)
        l2g2 = [m[2] for m in ms] .== Gene(0x2, 0x0)
        @test countnz(l1g1) / n == 1.0
        @test countnz(l2g2) / n == 1.0
        @test countnz(l1g1 & l2g2) / n == 1.0
    end

    ms = [meiosis(1, mito, 0.0, rec, Female, ch) for _ = 1:n]
    l1g1 = [m[1] for m in ms] .== Gene(0x1, 0x0)
    l2g2 = [m[2] for m in ms] .== Gene(0x2, 0x0)
    @test countnz(l1g1) / n == 1.0
    @test countnz(l2g2) / n == 1.0
    @test countnz(l1g1 & l2g2) / n == 1.0
end

ms = [meiosis(1, mito, 0.5, rec, Female, Male) for _ = 1:n]
l1g1 = [m[1] for m in ms] .== Gene(0x1, 0x0)
l2g2 = [m[2] for m in ms] .== Gene(0x2, 0x0)
@test countnz(l1g1) / n == 1.0
@test countnz(l2g2) / n == 1.0
@test countnz(l1g1 & l2g2) / n == 1.0

ms = [meiosis(1, mito, 0.5, rec, Female, Female) for _ = 1:n]
l1g1 = [m[1] for m in ms] .== Gene(0x1, 0x0)
l2g2 = [m[2] for m in ms] .== Gene(0x2, 0x0)
@test_approx_eq_eps countnz(l1g1) / n 0.5 tol
@test_approx_eq_eps countnz(l2g2) / n 0.5 tol
@test_approx_eq_eps countnz(l1g1 & l2g2) / n 0.25 tol

@test fuse(a.loci1, a.loci2, typeof(a)) == a
@test fuse(x.loci1, x.loci2, typeof(x)) == x
@test fuse(y2.loci1, y.loci1, typeof(y)) == y
@test fuse(mito.loci1, mito2.loci1, typeof(mito)) == mito

for par = (Female, Male)
    ret = migrate!(1, makeancestor(Autosome{1}), rec, par)
    @test getid(ret.loci1[1], Migration) != 0x0
    @test getid(ret.loci2[1], Migration) != 0x0
end

ret = migrate!(1, makeancestor(XChromosome{1}), rec, Female)
@test getid(ret.loci1[1], Migration) != 0x0
@test getid(ret.loci2[1], Migration) != 0x0
ret = migrate!(1, makeancestor(XChromosome{1}), rec, Male)
@test getid(ret.loci1[1], Migration) != 0x0
@test getid(ret.loci2[1], Migration) == 0x0

ret = migrate!(1, makeancestor(YChromosome{1}), rec, Female)
@test getid(ret.loci1[1], Migration) == 0x0
ret = migrate!(1, makeancestor(YChromosome{1}), rec, Male)
@test getid(ret.loci1[1], Migration) != 0x0

ret = migrate!(1, makeancestor(Mitochondrion{1}), rec, Female)
@test getid(ret.loci1[1], Migration) != 0x0
ret = migrate!(1, makeancestor(Mitochondrion{1}), rec, Male)
@test getid(ret.loci1[1], Migration) == 0x0

@test getgenotype(a) == [(Gene(0x1, 0x0), Gene(0x3, 0x0)), (Gene(0x2, 0x0), Gene(0x4, 0x0))]
@test getgenotype(x) == [(Gene(0x1, 0x0), Gene(0x3, 0x0)), (Gene(0x2, 0x0), Gene(0x4, 0x0))]
@test getgenotype(y) == [(Gene(0x1, 0x0),), (Gene(0x2, 0x0),)]
@test getgenotype(y2) == [(Gene(0x5, 0x0),), (Gene(0x6, 0x0),)]
@test getgenotype(mito) == [(Gene(0x1, 0x0),), (Gene(0x2, 0x0),)]
@test getgenotype(mito2) == [(Gene(0x3, 0x0),), (Gene(0x4, 0x0),)]
