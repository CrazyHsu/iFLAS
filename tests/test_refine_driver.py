"""End-to-end test of the refine LR driver on synthetic data.

Chain exercised:
  ground-truth isoforms -> cupcake-style collapsed GFF
    -> refine_isoforms (gtfToGenePred -> bed12+ -> strandAdjust -> readsAssign
       -> isoformGrouped.bed12+ + anno/novel split)
    -> assert all designed isoforms are recovered and each assigned to its gene.

ref_gpe / ref_bed are built from the synthetic genes.gtf (gtfToGenePred +
GenePredObj.toBed). TGS-only path (no NGS junctions). Requires bedtools
(pybedtools) + gtfToGenePred — both in the core `iflas` env.
"""

import os
import shutil
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import make_synthetic_data as syn  # noqa: E402

from iflas.cli import main as cli_main
from iflas.io.genepred import GenePredObj
from iflas.stages import refine as refine_stage
from iflas.stages.refine import refine_isoforms


pytestmark = pytest.mark.skipif(
    shutil.which("gtfToGenePred") is None or shutil.which("bedtools") is None,
    reason="needs gtfToGenePred + bedtools (core iflas env)")


def _write_collapsed_gff(path):
    with open(path, "w") as f:
        for gene, isos in syn.GENES.items():
            for name, chrom, strand, exons in isos:
                attrs = f'gene_id "{gene}"; transcript_id "{name}";'
                f.write(f"{chrom}\tPacBio\ttranscript\t{exons[0][0] + 1}\t{exons[-1][1]}\t.\t{strand}\t.\t{attrs}\n")
                for (st, en) in exons:
                    f.write(f"{chrom}\tPacBio\texon\t{st + 1}\t{en}\t.\t{strand}\t.\t{attrs}\n")


def _build_refs(gtf, ref_gpe, ref_bed):
    subprocess.run(["gtfToGenePred", "-genePredExt", gtf, ref_gpe], check=True)
    with open(ref_bed, "w") as f:
        for line in GenePredObj(ref_gpe, bincolumn=False).toBed(gene=True):
            f.write(str(line) + "\n")


def _rows(path):
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return [l.rstrip("\n").split("\t") for l in f if l.strip()]


def _write_wrapped_fasta(path, chrom, seq):
    with open(path, "w") as f:
        f.write(f">{chrom}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i + 80] + "\n")


def _bed12_junction(chrom, donor, acceptor, strand, score=20):
    chrom_start = donor - 1
    chrom_end = acceptor + 1
    return "\t".join(map(str, [
        chrom, chrom_start, chrom_end, f"junc_{donor}_{acceptor}", score, strand,
        chrom_start, chrom_end, "0,0,0", 2, "1,1", f"0,{acceptor - chrom_start}",
    ])) + "\n"


def _expected_isoform_count():
    return sum(len(isos) for isos in syn.GENES.values())


def _expected_gene_counts():
    return {gene: len(isos) for gene, isos in syn.GENES.items()}


def test_refine_recovers_isoforms_and_gene_assignment(tmp_path):
    d = str(tmp_path)
    syn.generate(d)
    gff = os.path.join(d, "collapsed.gff")
    _write_collapsed_gff(gff)
    ref_gpe = os.path.join(d, "ref.gpe")
    ref_bed = os.path.join(d, "ref.bed12+")
    _build_refs(os.path.join(d, "genes.gtf"), ref_gpe, ref_bed)

    out = refine_isoforms(gff, f"{d}/genome.fa", ref_gpe, ref_bed, d)

    grouped = _rows(out["isoformGrouped"])
    # all designed isoforms recovered
    assert len(grouped) == _expected_isoform_count()
    # each isoform assigned to its own gene (col4 = isoform name, col13 = gene)
    for row in grouped:
        assert row[12] == row[3].split(".")[0]
    # per-gene isoform counts match the design
    from collections import Counter
    counts = Counter(row[12] for row in grouped)
    assert counts == _expected_gene_counts()

    # exon structure preserved (geneB.iso1 = 3 exons 200,150,200)
    by_name = {row[3]: row for row in grouped}
    assert by_name["geneB.iso1"][10] == "200,150,200"
    assert by_name["geneB.iso2_SE"][9] == "2"          # SE isoform: 2 blocks
    assert by_name["geneD.iso2_IR"][9] == "1"          # IR isoform: single exon

    # anno/novel partition covers all isoforms
    anno = _rows(out["anno"])
    novel = _rows(out["novel"])
    assert len(anno) + len(novel) == _expected_isoform_count()
    # multi-exon isoforms whose intron chain matches the reference are annotated
    anno_names = {r[3] for r in anno}
    assert "geneB.iso1" in anno_names and "geneC.iso1" in anno_names


def test_refine_cli(tmp_path):
    d = str(tmp_path)
    syn.generate(d)
    gff = os.path.join(d, "collapsed.gff")
    _write_collapsed_gff(gff)
    ref_gpe = os.path.join(d, "ref.gpe")
    ref_bed = os.path.join(d, "ref.bed12+")
    _build_refs(os.path.join(d, "genes.gtf"), ref_gpe, ref_bed)

    rc = cli_main(["refine", "--collapsed-gff", gff, "--ref-genome", f"{d}/genome.fa",
                   "--ref-gpe", ref_gpe, "--ref-bed", ref_bed, "-o", d])
    assert rc == 0
    assert len(_rows(os.path.join(d, "refine", "isoformGrouped.bed12+"))) == _expected_isoform_count()


def test_refine_ngs_junctions_adjust_fuzzy_lr_junction(tmp_path):
    """NGS junction refinement replaces a supported fuzzy LR intron."""
    d = str(tmp_path)
    genome = list("A" * 1000)
    # original LR junction 200-500 and NGS-supported junction 205-495 are both
    # canonical; the replacement changes exon length by exactly 5%, so legacy
    # refineWithJunc accepts it.
    genome[200:202] = list("GT")
    genome[498:500] = list("AG")
    genome[205:207] = list("GT")
    genome[493:495] = list("AG")
    genome_fa = os.path.join(d, "genome.fa")
    _write_wrapped_fasta(genome_fa, "chr1", "".join(genome))

    gff = os.path.join(d, "collapsed.gff")
    with open(gff, "w") as f:
        attrs = 'gene_id "geneX"; transcript_id "pb_fuzzy";'
        f.write(f"chr1\tPacBio\ttranscript\t101\t600\t.\t+\t.\t{attrs}\n")
        f.write(f"chr1\tPacBio\texon\t101\t200\t.\t+\t.\t{attrs}\n")
        f.write(f"chr1\tPacBio\texon\t501\t600\t.\t+\t.\t{attrs}\n")

    ref_gpe = os.path.join(d, "ref.gpe")
    with open(ref_gpe, "w") as f:
        f.write("ref_tx\tchr1\t+\t100\t600\t100\t600\t2\t100,495,\t205,600,\t0\tgeneX\tcmpl\tcmpl\t0,0,\n")
    ref_bed = os.path.join(d, "ref.bed12+")
    with open(ref_bed, "w") as f:
        f.write("chr1\t100\t600\tref_tx\t0\t+\t100\t600\t0,0,0\t2\t105,105\t0,395\tgeneX\n")

    junctions = os.path.join(d, "junctions.bed")
    with open(junctions, "w") as f:
        f.write(_bed12_junction("chr1", 205, 495, "+", score=20))

    out = refine_isoforms(gff, genome_fa, ref_gpe, ref_bed, d, ngs_junctions=junctions, adjust=False)

    grouped = _rows(out["isoformGrouped"])
    assert len(grouped) == 1
    assert grouped[0][3] == "pb_fuzzy"
    assert grouped[0][9] == "2"
    assert grouped[0][10] == "105,105"
    assert grouped[0][11] == "0,395"
    assert grouped[0][12] == "geneX"
    assert _rows(out["anno"])[0][3] == "pb_fuzzy"

    adjusted = _rows(os.path.join(d, "refine", "tofu.juncAdjusted.bed12+"))
    assert adjusted[0][10] == "105,105"
    assert adjusted[0][11] == "0,395"


def test_refine_isoforms_uses_runner_for_gtf_to_genepred(tmp_path, monkeypatch):
    """The stage driver should route gtfToGenePred through ToolRunner."""
    d = str(tmp_path)
    gff = os.path.join(d, "collapsed.gff")
    genome = os.path.join(d, "genome.fa")
    ref_gpe = os.path.join(d, "ref.gpe")
    ref_bed = os.path.join(d, "ref.bed12+")
    for path in (gff, genome, ref_gpe, ref_bed):
        with open(path, "w") as f:
            f.write("")

    class RecordingRunner:
        def __init__(self):
            self.calls = []

        def run(self, cmd, **kwargs):
            self.calls.append((cmd, kwargs, os.getcwd()))
            with open("tofu.collapsed.gpe", "w") as f:
                f.write("fake\n")

    class DummyGpe:
        def __init__(self, path, bincolumn=False):
            self.path = path
            self.bincolumn = bincolumn

        def toBed(self, gene=True):
            assert self.path == "tofu.collapsed.gpe"
            assert gene is True
            return ["chr1\t0\t10\tiso1\t0\t+\t0\t10\t0,0,0\t1\t10\t0\tgene1"]

    def fake_reads_assign(ref_bed_arg, assign_input, readsColNum=13, outPrefix="tofu.collapsed.assigned", group=True):
        assert ref_bed_arg == os.path.abspath(ref_bed)
        assert assign_input == "tofu.collapsed.bed12+"
        with open("tofu.collapsed.assigned.unambi.bed12+", "w") as f:
            f.write("chr1\t0\t10\tiso1\t0\t+\t0\t10\t0,0,0\t1\t10\t0\tgroup1\tgene1\n")

    def fake_identify(ref_gpe_arg, bed_file, anno, novel):
        assert ref_gpe_arg == os.path.abspath(ref_gpe)
        assert bed_file == "isoformGrouped.bed12+"
        open(anno, "w").write("anno\n")
        open(novel, "w").write("")

    monkeypatch.setattr(refine_stage, "GenePredObj", DummyGpe)
    monkeypatch.setattr(refine_stage, "readsAssign", fake_reads_assign)
    monkeypatch.setattr(refine_stage, "identifyNovelIsoformsByJunctions", fake_identify)

    runner = RecordingRunner()
    out = refine_stage.refine_isoforms(gff, genome, ref_gpe, ref_bed, d, adjust=False, runner=runner)

    assert runner.calls == [
        (["gtfToGenePred", os.path.abspath(gff), "tofu.collapsed.gpe", "-genePredExt"],
         {}, os.path.join(d, "refine"))
    ]
    assert _rows(out["isoformGrouped"])[0][12] == "gene1"


def _cupcake_available():
    if shutil.which("conda") is None:
        return False
    r = subprocess.run(["conda", "run", "-n", "iflas-cupcake",
                        "collapse_isoforms_by_sam.py", "--help"], capture_output=True)
    return r.returncode == 0


@pytest.mark.skipif(not _cupcake_available(),
                    reason="iflas-cupcake env not available (cdna_cupcake)")
def test_real_collapse_to_refine_chain(tmp_path):
    """Full real chain in the sample layout: minimap2 -> junction filter (sorted BAM, addCVandID bed,
    processed.fa) -> cupcake collapse -> refine_isoforms incl. the reads-level assignment."""
    from iflas.stages.align import flnc_to_fasta, mapping_filter_and_add_tags, minimap2_tgs
    from iflas.stages.collapse import collapse_isoforms

    d = str(tmp_path)
    syn.generate(d)
    ref_gpe = os.path.join(d, "ref.gpe")
    ref_bed = os.path.join(d, "ref.bed12+")
    _build_refs(os.path.join(d, "genes.gtf"), ref_gpe, ref_bed)
    mapping = os.path.join(d, "mapping")
    sam = minimap2_tgs(f"{d}/genome.fa", f"{d}/flnc.fa", os.path.join(mapping, "flnc.mm2.sam"), platform="pacbio", threads=2)
    mapped = mapping_filter_and_add_tags(sam_file=sam, ref_bed_file=ref_bed, threads=2, work_dir=mapping)
    flnc_to_fasta(f"{d}/flnc.fa", os.path.join(mapping, "flnc.processed.fa"))
    os.makedirs(os.path.join(d, "collapse"))
    collapse_isoforms(f"{d}/flnc.fa", mapped["bam"], os.path.join(d, "collapse", "tofu"), fl_coverage=2)
    gff = os.path.join(d, "collapse", "tofu.collapsed.good.gff")

    out = refine_isoforms(gff, f"{d}/genome.fa", ref_gpe, ref_bed, d)
    grouped = _rows(out["isoformGrouped"])
    # cupcake yields distinct isoforms; all core truth genes remain assigned.
    from collections import Counter
    counts = Counter(row[12] for row in grouped)
    assert {"geneA", "geneB", "geneC", "geneD", "geneE"} <= set(counts)
    assert sum(counts.values()) == len(grouped) >= 7

    # reads level (legacy refineJunc "for reads"): every collapsed read assigned, exported as GTF
    used = set()
    for row in _rows(os.path.join(d, "collapse", "tofu.collapsed.group.txt")):
        used.update(row[1].split(","))
    reads = _rows(out["reads_assigned_unambi"])
    assert {row[3] for row in reads} == used and all(len(row) == 15 for row in reads)
    assert {row[14] for row in reads} <= set(counts)          # read genes == isoform genes
    gtf_tx = {l.split('transcript_id "')[1].split('"')[0] for l in open(out["reads_gtf"]) if "\ttranscript\t" in l}
    assert gtf_tx == used
