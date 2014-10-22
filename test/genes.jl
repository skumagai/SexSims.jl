getnewrec() = GeneStateRecorder(4)

rec = getnewrec()

@test rec.next == 0
@test rec.chunk == 4

gs = [mutate!(1, Gene(0, 0), rec) for _ = 1:4]
@test rec.next == 4
@test rec.chunk == 4

for (i, g) = enumerate(gs)
    @test id(g) == i
    @test mig(g) == 0
end

for (i, c) = enumerate(rec)
    @test c.gen == 1
    @test c.from == 0
    @test c.to == i
    @test c.kind == SexSims.MUTATION
end

mutate!(1, Gene(0, 0), rec)
@test rec.next == 5
@test rec.chunk == 4

rec = getnewrec()

gs = [migrate!(1, Gene(0, 0), rec) for _ = 1:4]
@test rec.next == 4
@test rec.chunk == 4

for (i, g) = enumerate(gs)
    @test id(g) == 0
    @test mig(g) == i
end

migrate!(1, Gene(0, 0), rec)
@test rec.next == 5
@test rec.chunk == 4

for (i, c) = enumerate(rec)
    @test c.gen == 1
    @test c.from == 0
    @test c.to == i
    @test c.kind == SexSims.MIGRATION
end
