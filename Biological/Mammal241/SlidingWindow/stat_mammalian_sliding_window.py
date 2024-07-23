import sys
import math
from scipy.stats import binom
import scipy.stats as stats
import numpy as np

def z_to_log_p(z, two_sided=True):
    # For large z, use log10(pdf(z)) - log10(z)
    if np.abs(z) > 12:  
        # Compute log10(pdf(x)) manually
        log_p = np.log10(1 / np.sqrt(2 * np.pi)) - (np.abs(z) ** 2 / 2) * np.log10(np.exp(1)) - np.log10(np.abs(z))
        if two_sided:
            # Multiply by 2 for two-sided p-value
            log_p += np.log10(2)
    else:
        if two_sided:
            # Multiply by 2 for two-sided p-value
            p = 2 * stats.norm.sf(np.abs(z))
        else:
            p = stats.norm.sf(z)
        # Use np.log10 to compute logarithm, adding a small value to avoid log10(0)
        log_p = np.log10(p + np.finfo(float).eps)
    return log_p

def median(arr):
    arr = sorted(arr)
    L = len(arr)
    if L % 2 == 0:
        return (arr[int(L / 2 - 1)] + arr[int(L / 2)]) / 2
    else:
        return arr[int(L / 2)]
        
def mean(arr):
    return sum(arr) / len(arr)

chrs = [str(i) for i in range(1, 23)]
nb = 1 + 237
BinSizes = [1, 5, 20, 100, 500]
cnt = [0 for i in range(len(BinSizes))]

QlowArr = [0 for i in range(nb)]

with open("stat_test_outlier.tsv", "w") as f:
    f.write("\t".join(["Branch", "Chr", "Pos", "Len", "Category", "Value", "relativeQ"]) + "\n")
    
    for branch in range(1, nb):
        A = [[] for i in range(len(BinSizes))]
        B = [[] for i in range(len(BinSizes))]
        C = [[] for i in range(len(BinSizes))]
        Q = [[] for i in range(len(BinSizes))]
        P = [[] for i in range(len(BinSizes))]
        a = {}
        b = {}
        c = {}
        q = {}
        p = {}
        n = {}
        for chromosome in chrs:
            headline = True
            a[chromosome] = []
            b[chromosome] = []
            c[chromosome] = []
            q[chromosome] = []
            p[chromosome] = []
            for line in open("sliding/" + str(branch) + "/chr" + chromosome + ".filter10k.tsv"):
                line = line.split()
                if headline:
                    headline = False
                    continue
                a[chromosome].append(float(line[1]))
                b[chromosome].append(float(line[2]))
                c[chromosome].append(float(line[3]))
                q[chromosome].append(float(line[4]))
                p[chromosome].append((chromosome, line[0]))
        sortedQ = sorted([e for chromosome in chrs for e in q[chromosome] if e > 0])
        Q1 = sortedQ[int(len(sortedQ) * 0.25 - 1e-6)]
        Q2 = sortedQ[int(len(sortedQ) * 0.5 - 1e-6)]
        Q3 = sortedQ[int(len(sortedQ) * 0.75 - 1e-6)]
        IQR = Q3 / Q1
        Qlow = Q1 / (IQR ** 1.5)
        QlowArr[branch] = Qlow
        for chromosome in chrs:
            n[chromosome] = [1 if q[chromosome][i] > Qlow else 0 for i in range(len(q[chromosome]))]
            a[chromosome] = [a[chromosome][i] if q[chromosome][i] > Qlow else 0 for i in range(len(q[chromosome]))]
            b[chromosome] = [b[chromosome][i] if q[chromosome][i] > Qlow else 0 for i in range(len(q[chromosome]))]
            c[chromosome] = [c[chromosome][i] if q[chromosome][i] > Qlow else 0 for i in range(len(q[chromosome]))]
            q[chromosome] = [q[chromosome][i] if q[chromosome][i] > Qlow else 0 for i in range(len(q[chromosome]))]
            for j in range(len(BinSizes)):
                bs = BinSizes[j]
                for i in range(int(len(q[chromosome]) / bs)):
                    if sum(n[chromosome][i*bs:(i+1)*bs]) >= bs * 0.5:
                        A[j].append(sum(a[chromosome][i*bs:(i+1)*bs]))
                        B[j].append(sum(b[chromosome][i*bs:(i+1)*bs]))
                        C[j].append(sum(c[chromosome][i*bs:(i+1)*bs]))
                        Q[j].append(sum(q[chromosome][i*bs:(i+1)*bs]))
                        P[j].append(p[chromosome][i*bs])
        for j in range(len(BinSizes)):
            bs = BinSizes[j]
            L = len(P[j])
            cnt[j] += L
            print(branch, bs, IQR, L)
            a = [A[j][i] / Q[j][i] for i in range(L)]
            b = [B[j][i] / Q[j][i] for i in range(L)]
            c = [C[j][i] / Q[j][i] for i in range(L)]
            q = [Q[j][i] for i in range(L)]
            p = [P[j][i] for i in range(L)]
            t = [a[i] + b[i] + c[i] for i in range(L)]
            sortedT = sorted(t)
            if bs < 100:
                T95 = sortedT[int(L * 0.95 - 1e-6)]
                T99 = sortedT[int(L * 0.99 - 1e-6)]
                for i in range(L):
                    v1 = np.log10(5) * (a[i] - T99) / (T99 - T95) - np.log10(0.01)
                    v2 = np.log10(5) * (b[i] - T99) / (T99 - T95) - np.log10(0.01)
                    if v1 > 4:
                        f.write("\t".join([str(branch), p[i][0], p[i][1], str(bs), "A", str(v1), str(q[i] / Q2 / bs)]) + "\n")
                    if v2 > 4:
                        f.write("\t".join([str(branch), p[i][0], p[i][1], str(bs), "B", str(v2), str(q[i] / Q2 / bs)]) + "\n")
            else:
                midT = sortedT[int(len(t) * 0.5 - 1e-6)]
                madT = sorted([abs(e - midT) for e in t])[int(len(t) * 0.5 - 1e-6)]
                for i in range(L):
                    if a[i] > midT:
                        v1 = -z_to_log_p(0.6745 * (a[i] - midT) / madT, two_sided=False)
                        if v1 > 3.5:
                            f.write("\t".join([str(branch), p[i][0], p[i][1], str(bs), "A", str(v1), str(q[i] / Q2 / bs)]) + "\n")
                    if b[i] > midT:
                        v2 = -z_to_log_p(0.6745 * (b[i] - midT) / madT, two_sided=False)
                        if v2 > 3.5:
                            f.write("\t".join([str(branch), p[i][0], p[i][1], str(bs), "B", str(v2), str(q[i] / Q2 / bs)]) + "\n")

threshold = 4
print("cnt:", cnt)
print("threshold:", threshold)
#exit(0)

with open("sliding.tsv", "w") as f2:
    buffer = []
    fcnt = 0
    f1cnt = 0
    binCount = {}
    f2.write("\t".join(["Branch", "Chr", "Pos", "Topology", "Score", "Coverage"]) + "\n")
    for branch in range(1, nb):
        for chromosome in chrs + ["X"]:
            headline = True
            cnt = 0
            A = 0
            B = 0
            C = 0
            Q = 0
            for line in open("sliding/" + str(branch) + "/chr" + chromosome + ".filter10k.tsv"):
                line = line.split()
                if headline:
                    headline = False
                    #outgroup = line
                    continue
                pos = int(line[0])
                s1 = float(line[1])
                s2 = float(line[2])
                s3 = float(line[3])
                q = float(line[4])
                if q > QlowArr[branch]:
                    A += s1
                    B += s2
                    C += s3
                    Q += q
                    cnt += 1
                offset = 4990000
                if pos % (offset + 10000) == offset:
                    if cnt < 250:
                        f2.write("\t".join([str(branch), chromosome, str(pos - offset), "A", "NaN", "low"]) + "\n")
                        f2.write("\t".join([str(branch), chromosome, str(pos - offset), "B", "NaN", "low"]) + "\n")
                        f2.write("\t".join([str(branch), chromosome, str(pos - offset), "C", "NaN", "low"]) + "\n")
                    else:
                        f2.write("\t".join([str(branch), chromosome, str(pos - offset), "A", str(A / Q), "high"]) + "\n")
                        f2.write("\t".join([str(branch), chromosome, str(pos - offset), "B", str(B / Q), "high"]) + "\n")
                        f2.write("\t".join([str(branch), chromosome, str(pos - offset), "C", str(C / Q), "high"]) + "\n")
                    A = 0
                    B = 0
                    C = 0
                    Q = 0
                    cnt = 0
                