import msprime
import numpy as np
import sys

PI = [18/61, 13/61, 14/61, 16/61]
GTR = msprime.GTR(relative_rates=[0.1, 0.3, 0.06, 0.12, 0.32, 0.1], equilibrium_frequencies=PI)

#(((A:0.2,B:0.08):0.001,C:0.04):0.1,D:0.3);
#(((A:20,B:8):0.1,C:4):10,D:30);

def AB_CD():
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=1, default_sampling_time = 0.0+0)
    demography.add_population(name="B", initial_size=1, default_sampling_time = 0.0+12)
    demography.add_population(name="C", initial_size=1, default_sampling_time = 0.1+16)
    demography.add_population(name="D", initial_size=1, default_sampling_time = 0.1+0)
    demography.add_population(name="I", initial_size=1)
    demography.add_population(name="ABC", initial_size=1)
    demography.add_population(name="ABCD", initial_size=1)
    demography.add_population_split(time=20, derived=["A", "B"], ancestral="I")
    demography.add_population_split(time=20.1, derived=["I", "C"], ancestral="ABC")
    demography.add_population_split(time=30.1, derived=["ABC", "D"], ancestral="ABCD")
    return demography

def AC_BD():
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=1, default_sampling_time = 0.0+0)
    demography.add_population(name="B", initial_size=1, default_sampling_time = 0.1+12)
    demography.add_population(name="C", initial_size=1, default_sampling_time = 0.0+16)
    demography.add_population(name="D", initial_size=1, default_sampling_time = 0.1+0)
    demography.add_population(name="I", initial_size=1)
    demography.add_population(name="ABC", initial_size=1)
    demography.add_population(name="ABCD", initial_size=1)
    demography.add_population_split(time=20, derived=["A", "C"], ancestral="I")
    demography.add_population_split(time=20.1, derived=["I", "B"], ancestral="ABC")
    demography.add_population_split(time=30.1, derived=["ABC", "D"], ancestral="ABCD")
    return demography

def BC_AD():
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=1, default_sampling_time = 0.1+0)
    demography.add_population(name="B", initial_size=1, default_sampling_time = 0.0+12)
    demography.add_population(name="C", initial_size=1, default_sampling_time = 0.0+16)
    demography.add_population(name="D", initial_size=1, default_sampling_time = 0.1+0)
    demography.add_population(name="I", initial_size=1)
    demography.add_population(name="ABC", initial_size=1)
    demography.add_population(name="ABCD", initial_size=1)
    demography.add_population_split(time=20, derived=["B", "C"], ancestral="I")
    demography.add_population_split(time=20.1, derived=["I", "A"], ancestral="ABC")
    demography.add_population_split(time=30.1, derived=["ABC", "D"], ancestral="ABCD")
    return demography



def simulate(demography, seed = 1, seqlen = 1000000):
    ts = msprime.sim_ancestry(samples={"A": 1, "B": 1, "C": 1, "D": 1}, demography=demography, ploidy=1, recombination_rate=0.001, sequence_length=seqlen, random_seed=seed)
    mts = msprime.sim_mutations(ts, rate=0.01, model=GTR, random_seed=seed)
    lst = mts.as_fasta().split("\n")
    L = len(lst) // 4
    default = np.random.choice(["A", "C", "G", "T"], seqlen, p=PI)
    aln = ["".join(lst[i*L+1:(i+1)*L]) for i in range(4)]
    return ["".join([default[j] if aln[i][j] == "N" else aln[i][j] for j in range(seqlen)]) for i in range(4)]

rep = int(sys.argv[1])
alns = [simulate(AB_CD(), seed = rep * 3 * 20 + 3 * i + 1) for i in range(20)]
for i in [1, 2, 5, 10, 20]:
    with open("ABCD_" + str(i) + "_" + str(rep) + ".fa", "w") as f:
        concat = ["".join([alns[j][k] for j in range(i)]) for k in range(4)]
        f.write("\n".join([">A", concat[0], ">B", concat[1], ">C", concat[2], ">D", concat[3], ""]))
alns = [simulate(AC_BD(), seed = rep * 3 * 20 + 3 * i + 2) for i in range(20)]
for i in [1, 2, 5, 10, 20]:
    with open("ACBD_" + str(i) + "_" + str(rep) + ".fa", "w") as f:
        concat = ["".join([alns[j][k] for j in range(i)]) for k in range(4)]
        f.write("\n".join([">A", concat[0], ">B", concat[1], ">C", concat[2], ">D", concat[3], ""]))
alns = [simulate(BC_AD(), seed = rep * 3 * 20 + 3 * i + 3) for i in range(20)]
for i in [1, 2, 5, 10, 20]:
    with open("BCAD_" + str(i) + "_" + str(rep) + ".fa", "w") as f:
        concat = ["".join([alns[j][k] for j in range(i)]) for k in range(4)]
        f.write("\n".join([">A", concat[0], ">B", concat[1], ">C", concat[2], ">D", concat[3], ""]))