# BioE131-231
**Introduction to Computational Biology (Fall 2019)**
*University of California, Berkeley*

A series of labs that focuses on fundamental bioinformatic practices.

**Lab 2:** Generates a phylogenetic tree (Newich-formatted tree) using sequences provided by the University using MUSCLE and FastTree. Identifies the sequences in the dataset using BLAST. Calculates sequences statistics for each identified cluster in the phylogenetic tree.

**Lab 3:** Takes the text "Romeo and Juliet" and plots a histogram of the words in the play.

**Lab 4:** Constructs a SQLite Database. Focuses on three pathways: glycolysis, the citrate cycle, and the pentose phosphate pathway. Considers genes in these pathways from Drosophila, *E. coli*, and humans when constructing associative tables.

**Lab 5:** Replicates the results from the publication “Computational design and experimental validation of oligonucleotide-sensing allosteric ribozymes” by Penchovsky and Breaker. Analyzes the ribozymes published in the paper using RNAfold.

**Lab 6:** Uses Illumina MiSeq reads from a sequencing run of *Shewanella oneidensis* to generate a plot showing how many times each position in the reference genome is covered. First it filters the reads for contamination from human cells.

**Lab 7:** Generates random binary and FASTA files which are then compressed using gzip (similar to LZW) and bzip2 (uses Burrows-Wheeler). Generates a table comparing the original file size, compressed file size, and run time for each algorithm on each sequence.

**Lab 8/9:** Using SPAdes, assembles a bacterial genome *de novo* using a combination of long PacBio reads and short Illumina reads. Analyzes the bacterial genome assembled using a few different statistics, identifies the genome’s taxonomic origin, and then obtains two genome annotations via different pipelines.

**Lab 10:**  Implements a verison of Needleman-Wunsch (NW) algorithm.

**Final Project:** This project’s original focus was to identify, tally, and compare Cas9 homologs including ISC genes within the transposons of ​*Salmonella​* genomes. However, we discovered that *​Salmonella* commonly uses the Type 1 CRISPR-Cas system. The components of which include: Cas6, Cas7, Cas5, Cas8, Cas3, Cas1, Cas2, and Cas4. Since our and other group’s genomes did not contain the sequences that encode Cas9 proteins (and consequently, any Cas9 homologs), we chose to study the presence of Type 1 Cas3 protein in the ​*Salmonella​* genomes.
