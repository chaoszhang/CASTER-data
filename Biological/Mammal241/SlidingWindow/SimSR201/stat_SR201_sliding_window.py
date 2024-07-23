import sys
import math
from scipy.stats import binom
import scipy.stats as stats
import numpy as np

cnt = 0
binCount = {}

with open("stat_pq_SIM.tsv", "w") as f:
    f.write("rep\T95ranch\tp\tq\n")
    for rep in range(1, 51):
        for branch in range(1, 199):
            A = []
            B = []
            C = []
            Q = []
            headline = True
            for line in open("sliding/" + str(rep) + "/sliding/" + str(branch) + "/10k.tsv"):
                line = line.split()
                if headline:
                    headline = False
                    continue
                s1 = float(line[1])
                s2 = float(line[2])
                s3 = float(line[3])
                q = float(line[4])
                if q > 0:
                    A.append(s1)
                    B.append(s2)
                    C.append(s3)
                    Q.append(q)
            L = len(Q)
            A = [A[i] / Q[i] for i in range(L)]
            B = [B[i] / Q[i] for i in range(L)]
            C = [C[i] / Q[i] for i in range(L)]
            T = [A[i] + B[i] + C[i] for i in range(L)]
            sortedT = sorted(T)
            midT = sortedT[int(L / 2)]
            cnt += L
            
            T95 = sortedT[int(L * 0.95 - 1e-6)]
            T99 = sortedT[int(L * 0.99 - 1e-6)]
            
            print(rep, branch)
            for i in range(1000):
                f.write(str(rep) + "\t" + str(branch) + "\t" + str((0.2 ** ((sortedT[int(i * L / 1000)] - T99) / (T99 - T95))) * 0.01) + "\t" + str(1-i/1000) + "\n")
            for t in T:
                b = round((t - T99) / (T99 - T95) * np.log2(5) - np.log2(0.01), 2)
                if not b in binCount:
                    binCount[b] = 0
                binCount[b] += 1
                
with open("stat_test_all.tsv", "w") as f1:
    f1.write("\t".join(["Bin", "Count"]) + "\n")
    for b in binCount:
        f1.write("\t".join([str(b), str(binCount[b])]) + "\n")
    