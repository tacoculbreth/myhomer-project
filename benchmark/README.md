# Benchmark
Performanace of myhomer was compared to SEA, Simple Enrichment Analysis in MEME Suite, and to UCSD's HOMER.

## Benchmarking against SEA
Result of SEA analysis for transcription factor Klf4 is in `benchmark/Klf4/sea_results` directory. Results comprise of SEA HTML ouptut, TSV output, passing sequences, uploaded BED file and sequences from the BED file.

## Benchmarking against HOMER
HOMER results for known motifs enrichment is in `benchmark/Klf4/homer_results/knownResults.html`.

To compare results of myhomer, check the output in `output/Klf4/motifs_enrichment_results.txt` for output of motifs in logo representation and PWM format. 

## Runtime Comparison
- myhomer: 163.96 seconds
- SEA: 34.43 seconds
- HOMER: 52 min 55 seconds
