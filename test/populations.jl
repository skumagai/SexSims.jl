using SexSims

facts("Test construction of Population") do
    context("10 organisms with four distinct types of chromosomes") do
        pop = Population(10, Autosome, XChromosome, YChromosome, Mitochondrion)
        @fact length(pop) => 10
        @fact chromosomes(pop) => \
            [Autosome, XChromosome, YChromosome, Mitochondrion]
        @fact stepsize(pop) => 6
        for i = 1:10
            @fact offset(pop, i) => (i - 1) * 6 + 1
        end
    end
    context("10 organisms with only x- and y-chromosomes") do
        pop = Population(10, XChromosome, YChromosome)
        @fact length(pop) => 10
        @fact chromosomes(pop) => [XChromosome, YChromosome]
        @fact stepsize(pop) => 3
        for i = 1:10
            @fact offset(pop, i) => (i - 1) * 3 + 1
    end
end

facts("Test migrations") do
    context("swap two individuals") do
        pop = Population(10, Autosome, Autosome)
        for i = 1:40
            pop.data[i] = Gene(i, i)
        end
        swap!(pop, 1, 2)
        for i = 1:4
            @fact pop.data[i] = Gene(i + 4, i + 4)
            @fact pop.data[i+4] = Gene(i, i)
        end
    end
end




