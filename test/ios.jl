params = SexSims.readinput(joinpath(dirname(@__FILE__()), "testinput.json"))
@test params.n == ([100, 200], [150, 150])
@test params.mut == 0.01
@test params.fit == ([1.0 0.5; 0.6 1.0], [1.0 0.7; 0.8 1.0])
@test params.l == (([0.01, 0.02], [0.02, 0.03]), ([0.1, 0.2], [0.2, 0.3]))
@test params.v == ([0.2 * 0.2 / (0.2 * 0.2 + 1 - 0.2), 0.3 * 0.1 / (0.3 * 0.1 + 1 - 0.3)], [0.4 * 0.02 / (0.4 * 0.02 + 1 - 0.4), 0.2 * 0.01 / (0.2 * 0.01 + 1 - 0.2)])
