# Sex of organisms
immutable Sex
    data::UInt8
end
const Female = Sex(0)
const Male   = Sex(1)

# Type of chromosomes
immutable ChromosomeKind
    data::UInt8
end
const Autosome      = ChromosomeKind(0)
const XChromosome   = ChromosomeKind(1)
const YChromosome   = ChromosomeKind(2)
const Mitochondrion = ChromosomeKind(3)

const ANCESTRAL = 0x0
