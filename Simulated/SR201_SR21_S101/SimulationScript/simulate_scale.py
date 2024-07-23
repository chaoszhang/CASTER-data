from dendropy.simulate import treesim
import dendropy
import msprime
import numpy as np
import collections
import sys

GENE_ID = "" if len(sys.argv) < 4 else "_" + sys.argv[3]
GEN_GAMMA_SHAPE = 5
LN_POP_SD = 0.5
LN_POP_MEAN = np.log(1e6) - (LN_POP_SD**2)/2
SEQ_LEN = int(5e6) if len(sys.argv) < 3 else int(sys.argv[2]) #2e7
MUT_RATE = 5e-8
LN_MUT_GAMMA_SHAPE_MEAN = 1.5
LN_MUT_GAMMA_SHAPE_SD = 1
MUT_LOCAL_LEN_MEAN = 1e4
MUT_GAMMA_SHAPE_LOCAL_MODIFIER = 20
REC_RATE = MUT_RATE * 1
EQ_FREQ_DIRICHLET = [36, 26, 28, 32]
GTR_MTX_DIRICHLET = [5, 15, 3, 6, 16, 5]

ORIGINAL_FILE_NAME = "model.200.10000000.0.0000001.strees"
FILE_GENERATION_MODIFIER = 0.1
CONDITION_FILE_NAME = "ASTERISK_S200_G1e6_B1e-7"

tree_list = dendropy.TreeList.get(path=CONDITION_FILE_NAME + ".strees", schema="newick")
rep = 0
for t in tree_list:
    rep += 1
    if rep != int(sys.argv[1]):
        continue
     #t = t.extract_tree_with_taxa_labels(["T0", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10",  "T100", "T200"])
    
    samplingTime = {}
    samples = []
    print("REP: " + str(rep))
    height = t.max_distance_from_root()
    for node in t.leaf_nodes():
        st = height - node.distance_from_root()
        node.edge_length += st
        samples.append(msprime.SampleSet(1, ploidy=1, time=st, population=node.taxon.label))
    tree = t.as_string(schema="newick", suppress_rooting=True)
    print("ancestry")
    initial_size = collections.defaultdict(lambda: np.random.lognormal(LN_POP_MEAN, LN_POP_SD))
    demography = msprime.Demography.from_species_tree(tree, initial_size)
    ts = msprime.sim_ancestry(samples=samples, ploidy=1, demography=demography, recombination_rate=REC_RATE, sequence_length=SEQ_LEN, random_seed=(hash(GENE_ID) % 100000)) #, model="smc_prime")
    for MUT_RATE_LEVEL in [1]:
        print("mutation")
        MUT_GAMMA_SHAPE = np.random.lognormal(LN_MUT_GAMMA_SHAPE_MEAN, LN_MUT_GAMMA_SHAPE_SD)
        MUT_GAMMA_SHAPE_LOCAL = MUT_GAMMA_SHAPE * MUT_GAMMA_SHAPE_LOCAL_MODIFIER
        rate = MUT_RATE_LEVEL*MUT_RATE*np.random.gamma(MUT_GAMMA_SHAPE_LOCAL, 1/MUT_GAMMA_SHAPE_LOCAL, size=SEQ_LEN)
        lastp = 0
        while lastp < SEQ_LEN:
            curp = min(lastp + np.random.geometric(1/MUT_LOCAL_LEN_MEAN) + 1, SEQ_LEN)
            rate[lastp:curp] *= np.random.gamma(MUT_GAMMA_SHAPE, 1/MUT_GAMMA_SHAPE)
            lastp = curp
        ratemap = msprime.RateMap(position=range(SEQ_LEN+1), rate=rate)
        EQ_FREQ = np.random.dirichlet(EQ_FREQ_DIRICHLET)
        while True:
            try:
                model = msprime.GTR(relative_rates=np.random.dirichlet(GTR_MTX_DIRICHLET), equilibrium_frequencies=EQ_FREQ)
                mts = msprime.sim_mutations(ts, rate=ratemap, model=model)
                break
            except:
                continue
        print("fasta")
        defaultSeq = np.random.choice(["A","C","G","T"], SEQ_LEN, p=EQ_FREQ)
        mts.write_fasta("data/original_fasta/" + str(rep) + "_" + str(MUT_RATE_LEVEL) + "X" + GENE_ID)
        with open("data/original_fasta/" + str(rep) + "_" + str(MUT_RATE_LEVEL) + "X" + GENE_ID, "r") as f:
            with open("data/alignment/" + str(rep) + "_" + str(MUT_RATE_LEVEL) + "X" + GENE_ID, "w") as f2:
                for line in f:
                    if line[0] == ">":
                        f2.write(">" + mts.populations()[mts.nodes()[int(line[2:])].population].metadata["name"] + "\n")
                        pos = 0
                    else:
                        for c in line:
                            f2.write(defaultSeq[pos] if c == "N" else c)
                            if c != "\n":
                                pos += 1

