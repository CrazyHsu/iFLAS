# Changelog

## [2.0.0b1] - 2026-10-16

First release of iFLAS for Python 3. It has been tested on simulated data, including reads simulated from the maize
B73 v5 reference, but not yet on real sequencing data. The numbers in brackets refer to our internal bug list and are
also used in comments in the code.

### Installation, configuration and commands

- iFLAS needs Python 3.10 or newer. All dependencies are installed from four conda environment files
  (`environment.yml`, `environment-cupcake.yml`, `environment-pacbio.yml`, `environment-diffsplice.yml`); the files in
  `envs/lock/` give the exact versions we tested.
- `config.yaml` and `samples.tsv` replace the .cfg file. `iflas migrate-config` converts an iFLAS 1 .cfg file. The
  settings are checked before a run starts.
- `iflas init-ref` makes the reference files (genePred, BED12, chromosome sizes, minimap2 and HISAT2 indexes) from a
  genome FASTA and a GTF or GFF3 file. Compressed files are accepted.
- `iflas all` runs the whole pipeline or the steps given with `--stages`; `-j` sets how many samples run at the same
  time. When a run is resumed without `preproc`, the reads from the earlier `preproc` are used.
- The options of the commands have long names (for example `--junc-comb-support`). The iFLAS 1 spellings (for example
  `-jcs`) are still accepted where possible. (BUG-129)
- New columns in the sample table: `fast5_dir` (nanopore signal for nanopolish). Nanopore reads can also be given as a
  dorado BAM file; its tail estimates (`pt` tag) are then used by palen_as.
- New settings: `pu_filter.random_state`, `options.ont_trim_polya`, `options.asas_ngs_confirm` and `dirs.tmp_dir`.
- Relative `out_dir` and `tmp_dir` are taken relative to the configuration file. (BUG-199)
- Gene structure plots (visual_as) are drawn with pyGenomeTracks instead of SpliceGrapher. R scripts are run with
  Rscript instead of rpy2. PyVCF and PyPDF2 are replaced by pysam and pypdf.
- The external tools are updated to current versions; see docs/external_tools.md.
- Large genomes: BAM files get CSI indexes when a chromosome is longer than 536,870,911 bp, and genomes larger than
  4 Gb get a single-part minimap2 index. (BUG-189, BUG-190)
- Paths may contain spaces. (BUG-168)
- visual_as is much faster on long genes: a 751-kb gene took over 40 CPU minutes and now takes under a second.
  (BUG-197)
- The report draws the allele-specific AS and APA plots of a sample in one R session instead of one per event; on our
  test data the report took 7 minutes instead of 27.

### Changes that affect results

These changes give different output from iFLAS 1. We measured the effect of each one on test data before adopting
it. Where a setting restores the iFLAS 1 behaviour, it is given.

Reads and alignments:

- Nanopore direct RNA reads: the basecalled 3' poly(A) is cut before mapping. With the tail left on, collapse lost
  about a third of the reads. The raw reads are still used by nanopolish. (BUG-178; `options.ont_trim_polya: false`)
- Nanopore reads are aligned against the genome FASTA, so that the nanopore settings of minimap2 (`-k14`) take
  effect. With the prebuilt index used by iFLAS 1 they were ignored. (BUG-192)
- The read filter after mapping compares each intron with `minimap2.max_intron_length`. iFLAS 1 compared the sum of
  all introns of a read, which removed all full-length reads of long genes. (BUG-193)
- nanopolish uses the alignments of the raw reads. iFLAS 1 combined the raw reads with alignments of the
  fmlrc2-corrected reads, and part of the tails were called SUFFCLIP. (BUG-175)
- HISAT2 runs with `--reorder`, so NGS alignments, StringTie assemblies and junction numbers are the same in every run.
  (BUG-163)
- When the reads of one strain are merged, FASTA and FASTQ files are converted to one format first. Before, the reads
  of one of the formats were lost in collapse. (BUG-177)

Isoforms:

- Isoforms with the same intron chain on opposite strands are kept apart, and each is labelled by the annotation of
  its own strand. (BUG-196)
- refine: splice-site dinucleotides are compared regardless of case, so soft-masked genomes are handled correctly
  (BUG-133); the tolerance for matching reference introns is the same on both sides (BUG-134); junction chains are
  matched as whole junctions (BUG-135, BUG-182); new gene clusters are formed per chromosome and strand (BUG-136);
  reads are assigned to genes on the same strand (BUG-137); the two thresholds of the strand correction are passed in
  the right order (BUG-142), and the strand evidence of single-exon reads (BUG-138) and of multi-exon reads without a
  clear splice motif (BUG-128) is used as intended; junction coordinates are sorted as numbers, not text (BUG-139);
  single-exon isoforms are labelled annotated or novel from their real overlap with the reference transcripts, where
  iFLAS 1 used a fixed test value (BUG-111).
- PU learning is seeded (`pu_filter.random_state`, default 0), so repeated runs give the same scores and the same
  high-quality isoforms. (BUG-141)
- The PU features on indels near junctions use the intended 0.5 ratio; iFLAS 1 used integer division. (BUG-140)
- `validIsoforms.lst` no longer contains a stray `pu_score` line, and low-quality annotated isoforms get a score.
  (BUG-179)
- When isoforms are grouped by junction chain (palen_as, asas), isoforms with one intron are no longer grouped with
  single-exon isoforms. (BUG-181)

AS events and functional analyses:

- asas and ase only use phased loci with at least one SNP that freebayes also finds in the NGS reads. IsoPhase
  sometimes reads alignment differences at splice sites as SNPs. (BUG-201; `options.asas_ngs_confirm: false` or
  `iflas asas --no-ngs-confirm`)
- Allele-specific expression works. In iFLAS 1 it never finished. (BUG-151) Allele counts from single-end NGS reads are
  no longer zero. (BUG-173)
- palen_as: in a merged strain, nanopore libraries with fast5 files keep their tail lengths when another library of
  the strain has a dorado tail table. (BUG-200)
- find_pa: reads ending at the right edge of a poly(A) cluster are counted in the cluster. (BUG-116)
- characterize: event lists are appended instead of overwritten (BUG-117), and the A3SS and A5SS splice-site tables
  are no longer swapped (BUG-118).
- diff_as: single-exon isoforms are no longer merged into one record per chromosome. (BUG-150)
- go: with several gene sets, the result table contains every set instead of only the last one (BUG-172); the plot is
  written to `<name>.goEnrichResults.pdf` (BUG-170).
- report: allele-specific AS plots of different gene clusters with the same reference gene no longer overwrite each
  other. (BUG-148)
- PacBio directories without a `subreadset.xml` file: `*.scraps.bam` is no longer used as subreads. (BUG-185)

### Other fixes

About sixty other problems, mostly crashes, are fixed. Among those that also occur in iFLAS 1: find_pa when no poly(A)
cluster passes the read threshold (BUG-153); characterize with splice sites near chromosome ends (BUG-154); asas and
ase on loci without phased SNPs, with more than two haplotypes, or with a haplotype without NGS reads (BUG-156,
BUG-161, BUG-180); palen_as with missing isoforms (BUG-157); merged strains in which some rows have no NGS reads
(BUG-176); fastp reports of several replicates overwriting each other (BUG-188); several report plots with current
R packages (BUG-164 to BUG-166, BUG-171); a clear error message for samples without NGS reads in pu_filter (BUG-191)
and for fast5 input without Guppy (BUG-187).

### Known limitations

- Not yet run on real sequencing data.
- Basecalling fast5 files needs Guppy (`guppy_basecaller`), which is no longer developed. nanopolish reads fast5 only;
  POD5 files have to be converted.
- Report plots that read alignments with Rsamtools are skipped for regions beyond 536.9 Mb.
- See section 8 of README.md for notes on the methods.
