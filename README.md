# Supplementary data of panModule

The panModule method is available in the PPanGGOLiN software suite: https://github.com/labgem/PPanGGOLiN

The article is currently unavailable.

## Benchmark

panModule predictions are benchmarked on a curated dataset of modules of _E. coli_ .

In the benchmark directory, you can find:
- 12_coli_IDs.tsv, a tsv file of the 12 _Escherichia coli_ strains with their assembly ids and accession numbers that are in the reference dataset
- EcoliComplet_genomes.list, a list of all 1671 genomes used for the EcoliComplet dataset.
- EcoliContigs_genomes.list, a list of all 1671 genomes used for the EcoliContigs dataset
- EcoliMAGs_benchmark_metrics.tsv, the results of the benchmark conducted with the EcoliMAGs dataset over all parameters.
- EcoliRef_benchmark_metrics.tsv, the results of the benchmark conducted with the EcoliScope dataset over all parameters.
- compute_modbench.py, the python script which computed the benchmark
- genomic_islands.tsv, the reference genomic islands in which the modules are found
- reference_modules.tsv, the reference modules


## _Klebsiella pneumoniae_ 1084

Analysis of _Klebsiella pneumoniae_ 1084 with panModule compared to [this study](
https://doi.org/10.1371/journal.pone.0096292).

the directory includes:
- Klebsiella_modules.tsv, a file listing the Klebsiella modules with a generic functional annotation and their positions on the genome.
- spot_19.html, an html file containing a dynamic figure of the spot_19, which contains the KPHPI208 region.
- spot_19_identical_rgps.tsv, a tsv file listing the identical RGPs of spot_19.
- Klebsiella_pneumoniae_genomes.list, a list of the 566 genomes used for the pangenome of Klebsiella pneumoniae.
