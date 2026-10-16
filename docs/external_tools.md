# External tools

This page lists the programs iFLAS calls, the versions we tested (the lock files in `envs/lock/`), and how the
command lines differ from those used by iFLAS 1.

## Tested versions

| Environment | Tools |
| :-- | :-- |
| iflas (core) | Python 3.11.15, minimap2 2.31, HISAT2 2.2.2, samtools 1.23.1, bedtools 2.31.1, StringTie 3.0.3, regtools 1.0.0, seqkit 2.13.0, fastp 1.3.3, fmlrc2 0.1.8, ropebwt2 r187, gffread 0.12.9, UCSC gtfToGenePred and genePredToGtf 482, nanopolish 0.14.0, freebayes 1.3.10, bcftools 1.24, bamtools 2.5.3, salmon 2.8.0, OpenJDK 25.0.2, pyGenomeTracks 3.9, Perl 5.32.1, scikit-learn 1.9.0, pysam 0.24.0, pandas 2.3.3, Biopython 1.87 |
| iflas-cupcake | Python 3.10.20, cDNA_Cupcake 29.0.0 (collapse and IsoPhase), PyVCF3 1.0.4, networkx 3.4.2 |
| iflas-pacbio | Python 3.7.12, ccs 6.4.0, lima 26.2.1, isoseq 26.2.0, pbbam 2.4.0 (pbindex), pbcoretools 0.8.1 (dataset) |
| iflas-diffsplice | rMATS 4.4.0, R 4.5.3 with Bioconductor 3.22: clusterProfiler 4.18.4, enrichplot 1.30.5, Gviz 1.54.0, Rsamtools 2.26.0, GenomicFeatures 1.62.0, txdbmaker 1.6.2; ggplot2 4.0.3, VennDiagram 1.8.2, dplyr 1.2.1, tidyr 1.3.2, valr 0.6.3 |

## Changes in the command lines

PacBio:

- ccs: versions 4 and later take `--min-passes`, `--min-length` and `--min-rq`, and no longer polish. `--min-length`
  is the minimum length of the CCS read, so iFLAS passes `ccs.min_ccs_length` there. `min_subread_length` and
  `min_read_score` are only used with ccs 3. When the input BAM already contains CCS or HiFi reads
  (`READTYPE=CCS` in the read group), ccs is skipped.
- lima names its output after the primer pair (`fl.5p--3p.bam`); iFLAS looks for the file lima actually wrote.
- isoseq3 was renamed isoseq in version 4, and iFLAS calls `isoseq`, so isoseq 4 or later is needed.
  `isoseq refine --require-polya` works as before.
- pbindex (from pbbam) indexes subreads BAM files that have no `.pbi` file before `dataset create`.
- bax2bam is only needed for RS II data (`.bax.h5`) and is not installed by default.

Nanopore:

- Basecalling is done outside iFLAS. Guppy is no longer developed; we suggest dorado. Give the FASTQ file, or the
  unaligned BAM file from dorado, as `data_location`, and the fast5 directory as `fast5_dir`.
- `nanopolish index` and `nanopolish polya` are called as in iFLAS 1. nanopolish reads fast5 only, so POD5 files have
  to be converted to fast5 first. For RNA004 data, the tail estimates of dorado (`--estimate-poly-a`) can be used
  instead.
- minimap2 aligns nanopore reads against the genome FASTA instead of a prebuilt index, so that `-k14` takes effect.

NGS:

- regtools 1.0 requires `-s` with a strand name. iFLAS uses `regtools junctions extract -s XS`, which gives the same
  result as `-s 0` in older versions (strand from the XS tag written by HISAT2).
- HISAT2 runs with `--reorder`, so the alignments come out in the order of the input reads.

Other tools:

- cDNA_Cupcake 29 names the collapse output `*.collapsed.good.gff` and exits with an error when no isoform is
  filtered out; iFLAS checks the output file instead of the exit code. The IsoPhase scripts are run with the Python of
  the cupcake environment, which also gets PyVCF3 and networkx (the conda package does not install them).
- rMATS 4.1 and later require `--tmp`, and since 4.0.2 reads with a length other than `--readLength` are dropped;
  iFLAS passes `--tmp` and `--variable-read-length`.
- salmon 2 writes `quant.sf`.
- `seqkit grep -f` is given `--id-regexp '^([^/\s]+)'` so that read names with `/1` and `/2` suffixes are matched.
- R scripts are run with Rscript (iFLAS 1 used rpy2). `makeTxDbFromGFF` comes from txdbmaker, because GenomicFeatures
  no longer provides it, and `show_guide` in ggplot2 calls is replaced by `show.legend`.
- Gene structure plots are drawn with pyGenomeTracks instead of SpliceGrapher. The conda package of pyGenomeTracks
  does not install tqdm, so `environment.yml` lists it.
- Biopython changed `IsoelectricPoint` after version 1.76. iFLAS keeps the 1.76 calculation for the isoelectric point
  used by CPC2, so the PU features are the same as in iFLAS 1.

## Effect of the newer versions

We rebuilt an environment with the tool versions of iFLAS 1 and ran both versions on the same simulated data. With
the same NGS alignments, the long-read alignments, isoform structures and PU features were the same. The newer HISAT2,
StringTie and regtools found somewhat different NGS junctions and junction read counts, which changes the features
that come from the NGS data. To compare with results from iFLAS 1 value by
value, use the NGS alignments and junctions of the iFLAS 1 run.
