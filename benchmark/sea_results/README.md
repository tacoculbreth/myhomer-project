# Benchmarking results from SEA, Simple Enrichment Analysis 

## Prameters for SEA Data Submission Form
- Type of control sequences to use: **shuffled input sequences**
- Sequence alphabet: **DNA, RNA or Protein**
- Input Sequence: **BED file**
- Genome: **UCSC Mammal Genomes, Mouse, mm10(GRCm38)**
- Input Motif: **MOUSE (Mus musculus) DNA, HOCOMOCO Mouse (v11 FULL)**

## How to prepare BED files (.bed) from HOMER Peaks file (peaks.txt)
```
tail -n +40 peaks.txt | cut -f 2-4 > Klf4_peaks.bed
awk '{print "chr"$0}' Klf4_peaks.bed > Klf4_peaks_ver2.bed
```

Make sure that .bed file format looks like:
```
chr17	35916928	35917002
chr17	36973463	36973537
chr17	57060668	57060742
chr17	35504026	35504100
chr17	45459121	45459195
chr17	57265117	57265191
chr17	83808182	83808256
chr17	87973887	87973961
chr17	85163102	85163176
```
