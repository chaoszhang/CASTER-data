import sys

seqs = {}

for line in open(sys.argv[1]):
    if len(line.split()) == 0:
        continue
    line = line.split()[0]
    if line[0] == '>':
        name = line
        seqs[name] = ""
    else:
        seqs[name] += line

H = ">homo_sap"
G = ">gallu_gal"
if H in seqs and G in seqs:
    pos = [i for i in range(len(seqs[H])) if seqs[H][i] in "acgtACGT" and seqs[G][i] in "acgtACGT"]
    for name in seqs:
        seqs[name] = [seqs[name][i] for i in pos]
    L = len(seqs[H])
    for i in range(L//50):
        dHG = len([j for j in range(i*50, (i+1)*50) if seqs[H][j] != seqs[G][j]])
        for name in seqs:
            d0 = len([j for j in range(i*50, (i+1)*50) if seqs[name][j] in "acgtACGT" and seqs[H][j] != seqs[name][j]])
            d1 = len([j for j in range(i*50, (i+1)*50) if seqs[name][j] in "acgtACGT" and seqs[H][j] != seqs[G][j]])
            if dHG > 10 or d0 > d1:
                for j in range(i*50, (i+1)*50):
                    seqs[name][j] = "-"
    OUTGROUP = [">gallu_gal", ">melea_gal", ">pelod_sin", ">anoli_car", ">xenop_tro", ">latim_cha", ">danio_rer", ">gaste_acu"]
    pos = [i for i in range((L//50)*50) if seqs[H][i] in "acgtACGT"]
    if len(pos) > 0:
        with open(sys.argv[1] + ".mask", "w") as f:
            for name in seqs:
                if name not in OUTGROUP:
                    f.write(name + "\n")
                    f.write("".join([seqs[name][i] for i in pos]) + "\n")
