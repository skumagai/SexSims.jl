using SexSims: Gene, GeneStateRecorder, getid, Mutation, Migration, generation, parent, kind, migrate!, mutate!, countalong, distance, eventintervals, listevents

n = 1000000
tol = 5e-2
srand(0x1)

@test Gene() == Gene(0x0, 0x0)

getnewrec() = GeneStateRecorder(4)

rec = getnewrec()

@test rec.next == 0
@test rec.chunk == 4

gs = [mutate!(1, Gene(0, 0), rec) for _ = 1:4]
@test rec.next == 4
@test rec.chunk == 4

for (i, g) = enumerate(gs)
    @test getid(g, Mutation) == i
    @test getid(g, Migration) == 0
end

for (i, c) = enumerate(rec)
    @test generation(c) == 1
    @test parent(c) == 0
    @test kind(c) == "mutation"
end

mutate!(1, Gene(0, 0), rec)
@test rec.next == 5
@test rec.chunk == 4

rec = getnewrec()

gs = [migrate!(1, Gene(0, 0), rec) for _ = 1:4]
@test rec.next == 4
@test rec.chunk == 4

for (i, g) = enumerate(gs)
    @test getid(g, Mutation) == 0
    @test getid(g, Migration) == i
end

migrate!(1, Gene(0, 0), rec)
@test rec.next == 5
@test rec.chunk == 4

for (i, c) = enumerate(rec)
    @test generation(c) == 1
    @test parent(c) == 0
    @test kind(c) == "migration"
end

rec = GeneStateRecorder(2)
gene = Gene()
gene = mutate!(1, gene, rec)
gene = mutate!(2, gene, rec)
gene = mutate!(3, gene, rec)
gene = mutate!(4, gene, rec)
gene = migrate!(5, gene, rec)
gene = migrate!(6, gene, rec)
gene = mutate!(7, gene, rec)
gene = migrate!(8, gene, rec)
gene = migrate!(9, gene, rec)
@test countalong(rec, gene, Migration) == 4
@test countalong(rec, gene, Mutation) == 5
@test sort(eventintervals(rec, gene, Mutation)) == [1, 1, 1, 3]
@test sort(eventintervals(rec, gene, Migration)) == [1, 1, 2]

rec = GeneStateRecorder(2)
gene = Gene()
gene1 = mutate!(1, gene, rec)
gene2 = mutate!(1, gene, rec)
gene1 = mutate!(2, gene1, rec)
gene1 = mutate!(3, gene1, rec)
gene2 = mutate!(2, gene2, rec)
@test listevents(rec, gene1, Mutation) == [4, 3, 1]
@test listevents(rec, gene2, Mutation) == [5, 2]
@test isnull(distance(rec, gene1, gene2, Mutation)) == true
gene1 = mutate!(1, gene, rec)
gene2 = mutate!(2, gene1, rec)
gene1 = mutate!(2, gene1, rec)
gene1 = mutate!(3, gene1, rec)
gene2 = mutate!(3, gene2, rec)
@test isnull(distance(rec, gene1, gene2, Mutation)) == false
@test get(distance(rec, gene1, gene2, Mutation)) == 5

rec = GeneStateRecorder(n)
gs = [getid(mutate!(1, Gene(), 0.5, rec), Mutation) for _ = 1:n]
@test_approx_eq_eps 0.5 (countnz(gs .== 0x0) / n) tol
gs = [getid(mutate!(1, Gene(), 0.75, rec), Mutation) for _ = 1:n]
@test_approx_eq_eps 0.25 (countnz(gs .== 0x0) / n) tol

