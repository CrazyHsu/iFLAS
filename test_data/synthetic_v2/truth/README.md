# synthetic_v2 truth tables (0-based, half-open coordinates unless stated)

| file | one row per | key columns |
|---|---|---|
| genes.tsv | gene | span, strand, category, softmasked, motif |
| isoforms.tsv | isoform (full set incl. novel) | exons, junc_chain, n_flnc, annotated (=in genes.gtf), junction_filter (legacy mapping filterByJunc outcome: pass / fail:chain-support / fail:consensus), op, wobble_shift, palen_mean |
| as_events.tsv | designed AS event | type, coords (find_as output convention), inc/exc isoforms, lr_psi, expected_detection in raw LR output (yes / no:junc-filter / no:R-2 / no:offset / no:MXE-unsupported / no:min-reads / no:aligner-microexon / both:A5SS+A3SS), expected_confident after the NGS confident filter (adds no:ngs-unsupported for alignment-wobble isoforms) |
| isoform_support.tsv | isoform | FLNC count and NGS fragment count per library |
| junction_support.tsv | designed junction | carrier isoforms, annotated, FLNC reads, designed NGS reads per library and in total (>= 5 bp anchors, as regtools -a 5 counts; the committed config merges all four libraries into one hisat2 run), ngs_detectable (0 when a flanking exon is shorter than the aligner can place, R-11), consensus for the legacy junction filter (annotated, or >= 3 FLNC reads, or > 10 NGS reads in total) |
| read_truth.tsv | FLNC read | isoform, haplotype, PA site, poly(A) length + nanopolish qc, jitter/truncation, wobble, errors, primers + zmw + `flnc_recovered` (what lima/isoseq refine keep from ccs.bam) |
| diff_psi.tsv | NGS differential event | designed PSI per condition, expected_significant |
| ase.tsv | gene with SNPs | designed FL reads per haplotype over all isoforms, hapA share, expected_imbalanced (|share - 0.5| >= 0.25: allele-specific expression) |
| allele_bias.tsv | ASAS event | P(hapA | inclusion), designed FL reads per haplotype of the inclusion/exclusion isoform (exact quota), expected_significant under the published criterion (>= 5 / <= 3 reads per haplotype per isoform, opposite majority haplotypes, chi-square p <= 0.001) |
| pa_sites.tsv | PA site | genomic site, fraction, designed reads |
| palen.tsv | AS event | designed poly(A) means of inc/exc isoforms, PASS reads per side, expected_significant (means differ by design), detectable_published (the published t-test + Kruskal <= 0.001 rule with >= 5 PASS reads per side can resolve it at the designed depth) |
| genome_features.tsv | interval | softmasked / noncanonical_splice / te_repeat |
| snps.vcf | heterozygous SNP | 1-based VCF, INFO GENE |
| genes.full.gtf | all isoforms | truth annotation (genes.gtf is the *reference subset* the pipeline sees) |
