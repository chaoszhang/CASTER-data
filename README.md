# CASTER: Direct species tree inference from whole-genome alignments
Scripts for generating data and figures

## Proof
Mathmatica script to verify the correctness of Proposition 2

## Simulated
Simulated datasets

### Quartet
Scripts for Quartet dataset

- SimulationScript: MSPrime scripts to simulate quartet phylogenies

- BPP_ControlFiles: Templates for BPP control files

- IntrogressionSignal: Implementation for parsimony-based introgression signal detection for quartets (ABBA-BABA test)

- FigureScripts: Script and data for generating Figure 2 and related figures

### SR201_SR21_S101
Scripts for SR201 dataset, SR21 dataset, and S101 datset

- SimulationScript: MSPrime scripts to simulate SR201 and SR21 phylogenies

- FigureScripts: Script and data for generating Figure 3 and related figures

## Biological
Scripts used in biological datasets

### Mammal241
Foley et al. placental mammal dataset

- Phylogeny: whole-genome phylogenies and per-chromosome phylogenies by CASTER-site

- SlidingWindow: Script and data for generating Figure 4 and related figures for sliding window analysis

### Bird363
Stiller et al. bird dataset

- Phylogeny: Input file list and output phylogenies of CASTER-site

- SlidingWindow: Script and data for generating Figure 5 and related figures

### Others
Relavent scripts used in other biological datasets

# Commands

CASTER (site and pair):

```
$PROGRAM_NAME -C -t $THREADS $MSA
```

wASTRAL (gene trees):

```
raxml-ng --msa $MSA --model GTR+G --threads 1
iqtree2 -s $MSA -te $MSA.raxml.bestTree -m GTR+G4 -abayes -nt 1
```

wASTRAL (species tree):

```
astral-hybrid -B -C -t $THREADS $GENE_TREES
```

ASTRAL (species trees):

```
astral -C -t $THREADS $GENE_TREES
```

SVDQuartets (haploid, single individual):

```
svdquartets nthreads=$THREADS evalQuartets=all seed=5000
```

SVDQuartets (unphased diploid):

```
svdquartets nthreads=$THREADS evalQuartets=random nquartets=105597360 taxpartition=species seed=5000
```

SVDQuartets (multiple individuals):

```
svdquartets nthreads=$THREADS evalQuartets=all taxpartition=species seed=5000
```

RAxML-ng (haploid, single individual):

```
raxml-ng --tree pars{1} --msa $MSA --model GTR+G --threads $THREADS
```

RAxML-ng (unphased diploid):

```
raxml-ng --tree pars{1} --msa $MSA --model GTGTR4+G --threads $THREADS
```

RAxML-ng (multiple individuals):

```
raxml-ng --tree pars{1} --msa $MSA --model GTR+G --threads $THREADS --tree-constraint $SPECIES_DELIMITATION
```

Placental mammal whole genome alignment:

```
caster-site -C -t 64 -o chr1.nw chr1.fa
caster-site -r 0 -s 0 -t 128 -f list -g chr1.nw -o whole_genome.nw chr_list.txt
```
