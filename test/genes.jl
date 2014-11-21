n = 1000000
tol = 5e-2
srand(0x1)

getnewrec() = GeneStateRecorder(4)

rec = getnewrec()

@test rec.next == 0
@test rec.chunk == 4

gs = [mutate!(1, Gene(0, 0), rec) for _ = 1:4]
@test rec.next == 4
@test rec.chunk == 4

for (i, g) = enumerate(gs)
    @test SexSims.getid(g, SexSims.Mutation) == i
    @test SexSims.getid(g, SexSims.Migration) == 0
end

for (i, c) = enumerate(rec)
    @test SexSims.generation(c) == 1
    @test SexSims.parent(c) == 0
    @test SexSims.kind(c) == "mutation"
end

mutate!(1, Gene(0, 0), rec)
@test rec.next == 5
@test rec.chunk == 4

rec = getnewrec()

gs = [migrate!(1, Gene(0, 0), rec) for _ = 1:4]
@test rec.next == 4
@test rec.chunk == 4

for (i, g) = enumerate(gs)
    @test SexSims.getid(g, SexSims.Mutation) == 0
    @test SexSims.getid(g, SexSims.Migration) == i
end

migrate!(1, Gene(0, 0), rec)
@test rec.next == 5
@test rec.chunk == 4

for (i, c) = enumerate(rec)
    @test SexSims.generation(c) == 1
    @test SexSims.parent(c) == 0
    @test SexSims.kind(c) == "migration"
end

rec = GeneStateRecorder(2)
gene = Gene(0, 0)
gene = mutate!(1, gene, rec)
gene = mutate!(2, gene, rec)
gene = mutate!(3, gene, rec)
gene = mutate!(4, gene, rec)
gene = migrate!(5, gene, rec)
gene = migrate!(6, gene, rec)
gene = mutate!(7, gene, rec)
gene = migrate!(8, gene, rec)
gene = migrate!(9, gene, rec)
@test SexSims.countalong(rec, gene, SexSims.Migration) == 4
@test SexSims.countalong(rec, gene, SexSims.Mutation) == 5
@test sort(SexSims.eventintervals(rec, gene, SexSims.Mutation)) == [1, 1, 1, 3]
@test sort(SexSims.eventintervals(rec, gene, SexSims.Migration)) == [1, 1, 2]

rec = GeneStateRecorder(2)
gene = Gene(0, 0)
gene1 = mutate!(1, gene, rec)
gene2 = mutate!(1, gene, rec)
gene1 = mutate!(2, gene1, rec)
gene1 = mutate!(3, gene1, rec)
gene2 = mutate!(2, gene2, rec)
@test SexSims.listevents(rec, gene1, SexSims.Mutation) == [4, 3, 1]
@test SexSims.listevents(rec, gene2, SexSims.Mutation) == [5, 2]
@test isnull(SexSims.distance(rec, gene1, gene2, SexSims.Mutation)) == true
gene1 = mutate!(1, gene, rec)
gene2 = mutate!(2, gene1, rec)
gene1 = mutate!(2, gene1, rec)
gene1 = mutate!(3, gene1, rec)
gene2 = mutate!(3, gene2, rec)
@test isnull(SexSims.distance(rec, gene1, gene2, SexSims.Mutation)) == false
@test get(SexSims.distance(rec, gene1, gene2, SexSims.Mutation)) == 5

rec = GeneStateRecorder(n)
gs = [SexSims.getid(mutate!(1, Gene(0, 0), 0.5, rec), SexSims.Mutation) for _ = 1:n]
@test_approx_eq_eps 0.5 (countnz(gs .== 0x0) / n) tol
gs = [SexSims.getid(mutate!(1, Gene(0, 0), 0.75, rec), SexSims.Mutation) for _ = 1:n]
@test_approx_eq_eps 0.25 (countnz(gs .== 0x0) / n) tol

