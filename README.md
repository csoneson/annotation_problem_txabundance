## A junction coverage compatibility score to quantify the reliability of transcript abundance estimates and annotation catalogs

This repository contains all the code that was used to perform the evaluations in the manuscript defining the JCC score:

Soneson C, Love MI, Patro R, Hussain S, Malhotra D and Robinson MD: [A junction coverage compatibility score to quantify the reliability of transcript abundance estimates and annotation catalogs](https://www.biorxiv.org/content/early/2018/07/28/378539). bioRxiv (2018).

The main purpose of this repository is to provide a record of the analyses in our paper. All the results were obtained using the provided `Makefile`. Note, however, that reproducing the results requires installing all the software tools that were used, downloading the raw RNA-seq data as well as the annotation files, and subsequently setting the correct paths to these (typically done in the beginning of the master `Makefile` or any of the specific makefiles in the `makefiles` subdirectory).

### Notes
The `fasta_to_fastq.pl` script was downloaded from [https://code.google.com/archive/p/fasta-to-fastq/](https://code.google.com/archive/p/fasta-to-fastq/) on May 17, 2018.
