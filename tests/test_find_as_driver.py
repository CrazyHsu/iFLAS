"""End-to-end test of the find_as LR driver on synthetic data.

Chain exercised:
  hqIsoforms/hq.collapsed.bed12+ + collapse/tofu.collapsed.group.txt
    -> find_as_sample / find_as_lr
    -> LR/{IR,SE,A5SS,A3SS} files, asserted against the designed AS events.

Reference (ref_gpe) is built from the synthetic genes.gtf via gtfToGenePred, and
the read<->gene assignment is done by the real getAnnotation overlap logic, so this
covers the realistic find_all_as LR path (minus the deferred NGS junction filtering).

Requires bedtools (pybedtools) + gtfToGenePred — both in the core `iflas` env.
"""

import os
import shutil
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import make_synthetic_data as syn  # noqa: E402

from iflas.cli import main as cli_main
from iflas.stages.find_as import (
    apply_ngs_confident_filters,
    collapsed_gff_to_bed12plus,
    find_as_lr,
    find_as_sample,
)


pytestmark = pytest.mark.skipif(
    shutil.which("gtfToGenePred") is None or shutil.which("bedtools") is None or shutil.which("perl") is None,
    reason="needs gtfToGenePred + bedtools + perl (core iflas env)")


def _write_collapsed_gff(path):
    """Write a cupcake-style collapsed GFF straight from the ground-truth isoforms.

    transcript_id == isoform name so the assertions can name expected reads; the
    group file (see ``_write_group``) maps each isoform to fake CCS reads.
    """
    with open(path, "w") as f:
        for gene, isos in syn.GENES.items():
            for name, chrom, strand, exons in isos:
                tx_start = exons[0][0] + 1   # GFF is 1-based inclusive
                tx_end = exons[-1][1]
                attrs = f'gene_id "{gene}"; transcript_id "{name}";'
                f.write(f"{chrom}\tPacBio\ttranscript\t{tx_start}\t{tx_end}\t.\t{strand}\t.\t{attrs}\n")
                for (st, en) in exons:
                    f.write(f"{chrom}\tPacBio\texon\t{st + 1}\t{en}\t.\t{strand}\t.\t{attrs}\n")


def _write_group(path, copies_per_iso=5):
    with open(path, "w") as f:
        for gene, isos in syn.GENES.items():
            for name, chrom, strand, exons in isos:
                reads = ",".join(f"{name}/ccs/{k}" for k in range(copies_per_iso))
                f.write(f"{name}\t{reads}\n")


def _write_hq_bed12plus(path):
    """Write the post-pu_filter bed12+ contract consumed by find_as.

    Column 15 (``otherList[2]``) is the collapse group key. The surrounding
    payload columns mirror the wider hq.collapsed.bed12+ shape but are not read
    by find_as.
    """
    with open(path, "w") as f:
        for gene, isos in syn.GENES.items():
            for name, chrom, strand, exons in isos:
                chrom_start, chrom_end = exons[0][0], exons[-1][1]
                block_sizes = ",".join(str(en - st) for st, en in exons)
                block_starts = ",".join(str(st - chrom_start) for st, en in exons)
                fields = [
                    chrom, chrom_start, chrom_end, name, 0, strand,
                    chrom_start, chrom_end, "0,0,0", len(exons),
                    block_sizes, block_starts,
                    gene, gene, name, len(exons), "annotated", 5, 5,
                ]
                f.write("\t".join(map(str, fields)) + "\n")


def _write_find_as_inputs(base_dir):
    hq_dir = os.path.join(base_dir, "hqIsoforms")
    collapse_dir = os.path.join(base_dir, "collapse")
    os.makedirs(hq_dir, exist_ok=True)
    os.makedirs(collapse_dir, exist_ok=True)
    hq_bed = os.path.join(hq_dir, "hq.collapsed.bed12+")
    group = os.path.join(collapse_dir, "tofu.collapsed.group.txt")
    _write_hq_bed12plus(hq_bed)
    _write_group(group)
    return hq_bed, group


def _build_ref_gpe(gtf, out_gpe):
    subprocess.run(["gtfToGenePred", "-genePredExt", gtf, out_gpe], check=True)


def _rows(path):
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return [l.rstrip("\n").split("\t") for l in f if l.strip()]


def _write_junctions(path, junctions):
    """Write simple BED12 junctions whose Junction(j).jPos is start-end."""
    with open(path, "w") as f:
        for chrom, donor, acceptor, strand, score in junctions:
            chrom_start = donor - 1
            block_starts = f"0,{acceptor - chrom_start}"
            fields = [
                chrom, chrom_start, acceptor + 1, f"junc_{chrom}_{donor}_{acceptor}",
                score, strand, chrom_start, acceptor + 1, "0,0,0", 2,
                "1,1", block_starts,
            ]
            f.write("\t".join(map(str, fields)) + "\n")


def test_bridge_gff_to_bed12plus(tmp_path):
    d = str(tmp_path)
    syn.generate(d)
    gff = os.path.join(d, "collapsed.gff")
    grp = os.path.join(d, "collapsed.group.txt")
    _write_collapsed_gff(gff)
    _write_group(grp)
    bed = os.path.join(d, "iso.bed12+")
    collapsed_gff_to_bed12plus(gff, grp, bed)

    by_name = {r[3]: r for r in _rows(bed)}
    assert set(by_name) == {iso[0] for isos in syn.GENES.values() for iso in isos}
    # geneB.iso1: 3 exons, 0-based 2000..3200, group key carried in col15
    r = by_name["geneB.iso1"]
    assert r[0] == "chr1" and r[1] == "2000" and r[2] == "3200"
    assert r[5] == "+" and r[9] == "3"
    assert r[10] == "200,150,200"          # block sizes
    assert r[11] == "0,500,1000"           # block rel starts
    assert r[14] == "geneB.iso1"           # otherList[2] -> group key
    # single-exon IR isoform
    assert by_name["geneD.iso2_IR"][9] == "1"


def test_find_as_lr_detects_all_event_types(tmp_path):
    d = str(tmp_path)
    syn.generate(d)
    base = os.path.join(d, "sample")
    bed, grp = _write_find_as_inputs(base)

    ref_gpe = os.path.join(d, "ref.gpe")
    _build_ref_gpe(os.path.join(d, "genes.gtf"), ref_gpe)

    out = find_as_lr(bed, grp, ref_gpe, d)

    # SE: geneB skips exon 2500-2650
    se = _rows(out["SE"])
    se_ids = [r[3] for r in se]
    assert "geneB:2200@2500-2650@3000" in se_ids
    se_row = next(r for r in se if r[3] == "geneB:2200@2500-2650@3000")
    assert se_row[0] == "chr1" and se_row[1] == "2500" and se_row[2] == "2650"
    assert se_row[12] == "geneB"

    # IR: geneD retains intron 300-700
    ir = _rows(out["IR"])
    ir_ids = [r[3] for r in ir]
    assert "geneD:300-700" in ir_ids
    ir_row = next(r for r in ir if r[3] == "geneD:300-700")
    assert ir_row[0] == "chr2" and ir_row[6] == "geneD"

    # A5SS: geneC alt donor 4000/4050
    a5 = _rows(out["A5SS"])
    a5_ids = [r[3] for r in a5]
    assert "geneC:4000-4050" in a5_ids

    # A3SS: geneE alt acceptor 4950/5000
    a3 = _rows(out["A3SS"])
    a3_ids = [r[3] for r in a3]
    assert "geneE:4950-5000" in a3_ids


def test_find_as_lr_writes_confident_links_without_ngs(tmp_path):
    d = str(tmp_path)
    syn.generate(d)
    base = os.path.join(d, "sample")
    bed, grp = _write_find_as_inputs(base)

    ref_gpe = os.path.join(d, "ref.gpe")
    _build_ref_gpe(os.path.join(d, "genes.gtf"), ref_gpe)

    out = find_as_lr(bed, grp, ref_gpe, d)

    for as_type in ("IR", "SE", "A5SS", "A3SS"):
        assert os.path.islink(out["confident"][as_type])
        assert _rows(out["confident"][as_type]) == _rows(out[as_type])


def test_find_as_lr_filters_ir_confident_by_ngs_junctions(tmp_path):
    d = str(tmp_path)
    syn.generate(d)
    base = os.path.join(d, "sample")
    bed, grp = _write_find_as_inputs(base)

    ref_gpe = os.path.join(d, "ref.gpe")
    _build_ref_gpe(os.path.join(d, "genes.gtf"), ref_gpe)
    junctions = os.path.join(d, "junctions.bed")
    _write_junctions(junctions, [("chr2", 300, 700, "+", 12)])

    out = find_as_lr(bed, grp, ref_gpe, d, ngs_junctions=junctions)

    assert [r[3] for r in _rows(out["confident"]["IR"])] == ["geneD:300-700"]
    assert os.path.islink(out["confident"]["IR"])
    assert os.path.basename(os.readlink(out["confident"]["IR"])) == "IR.NGS.bed6+"


def test_apply_ngs_confident_filters_uses_runner(tmp_path, monkeypatch):
    from iflas.stages import scan_se_by_ngs

    d = str(tmp_path)
    lr_dir = os.path.join(d, "LR")
    os.makedirs(lr_dir)
    out_files = {
        "IR": os.path.join(lr_dir, "IR.bed6+"),
        "SE": os.path.join(lr_dir, "SE.bed12+"),
        "A5SS": os.path.join(lr_dir, "A5SS.bed6+"),
        "A3SS": os.path.join(lr_dir, "A3SS.bed6+"),
    }
    with open(out_files["IR"], "w") as f:
        f.write("chr1\t10\t20\tgeneX:10-20\t0\t+\tgeneX\tir\t1\tspliced\t1\n")
    with open(out_files["SE"], "w") as f:
        f.write("chr1\t30\t40\tgeneX:se\t0\t+\t30\t40\t0,0,0\t1\t10\t0\tgeneX\n")
    for as_type in ("A5SS", "A3SS"):
        open(out_files[as_type], "w").close()

    isoform_bed = os.path.join(d, "hq.bed12+")
    with open(isoform_bed, "w") as f:
        f.write("chr1\t0\t100\tiso1\t0\t+\t0\t100\t0,0,0\t2\t10,10\t0,90\tgeneX\tgeneX\tiso1\n")
    ref_gpe = os.path.join(d, "ref.gpe")
    open(ref_gpe, "w").close()
    junctions = os.path.join(d, "junctions.bed")
    _write_junctions(junctions, [("chr1", 10, 20, "+", 12)])

    def fake_scan(**kwargs):
        with open(kwargs["outFile"], "w") as f:
            f.write("chr1\t30\t40\tgeneX:se\n")

    def fake_assign(**kwargs):
        with open(kwargs["outFile"], "w") as f:
            f.write("chr1\t30\t40\tgeneX:se\n")
        open(kwargs["errorOutFile"], "w").close()

    monkeypatch.setattr(scan_se_by_ngs, "scanSeByNGS", fake_scan)
    monkeypatch.setattr(scan_se_by_ngs, "assignGeneNameToSE", fake_assign)

    class RecordingRunner:
        def __init__(self):
            self.calls = []

        def run(self, cmd, **kwargs):
            self.calls.append((cmd, kwargs))
            script = os.path.basename(cmd[1])
            stdout = kwargs.get("stdout")
            stderr = kwargs.get("stderr")
            if script == "juncAssign.pl":
                stdout.write("chr1\t9\t21\tj1\t12\t+\t9\t21\t0,0,0\t2\t1,1\t0,11\n")
            elif script in {"AnSSconfirmByJunc.pl", "getNgsA5SS.pl", "getNgsA3SS.pl"}:
                if kwargs.get("input") is not None:
                    assert "j1:+" in kwargs["input"]
                stdout.write("")
                if stderr is not None:
                    stderr.write("")
            else:
                raise AssertionError(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    runner = RecordingRunner()
    confident = apply_ngs_confident_filters(
        out_files, isoform_bed, ref_gpe, junctions, d, runner=runner
    )

    assert [os.path.basename(call[0][1]) for call in runner.calls] == [
        "juncAssign.pl",
        "AnSSconfirmByJunc.pl",
        "AnSSconfirmByJunc.pl",
        "getNgsA5SS.pl",
        "getNgsA3SS.pl",
    ]
    assert all("stdout" in kwargs for _, kwargs in runner.calls)
    assert os.path.islink(confident["IR"])


def test_find_as_sample_uses_hq_bed_contract(tmp_path):
    d = str(tmp_path)
    syn.generate(d)
    base = os.path.join(d, "sample")
    _write_find_as_inputs(base)
    ref_gpe = os.path.join(d, "ref.gpe")
    _build_ref_gpe(os.path.join(d, "genes.gtf"), ref_gpe)

    out = find_as_sample(base, ref_gpe)
    assert out["SE"] == os.path.join(base, "as_events", "ordinary_as", "LR", "SE.bed12+")
    assert "geneB:2200@2500-2650@3000" in [r[3] for r in _rows(out["SE"])]
    assert "geneD:300-700" in [r[3] for r in _rows(out["IR"])]


def test_find_as_cli_base_dir_mode(tmp_path):
    d = str(tmp_path)
    syn.generate(d)
    base = os.path.join(d, "sample")
    _write_find_as_inputs(base)
    ref_gpe = os.path.join(d, "ref.gpe")
    _build_ref_gpe(os.path.join(d, "genes.gtf"), ref_gpe)

    rc = cli_main(["find_as", "-b", base, "-r", ref_gpe])
    assert rc == 0
    lr_dir = os.path.join(base, "as_events", "ordinary_as", "LR")
    assert "geneB:2200@2500-2650@3000" in [r[3] for r in _rows(os.path.join(lr_dir, "SE.bed12+"))]
    assert "geneD:300-700" in [r[3] for r in _rows(os.path.join(lr_dir, "IR.bed6+"))]


def _cupcake_available():
    if shutil.which("conda") is None:
        return False
    r = subprocess.run(["conda", "run", "-n", "iflas-cupcake",
                        "collapse_isoforms_by_sam.py", "--help"], capture_output=True)
    return r.returncode == 0


@pytest.mark.skipif(not _cupcake_available(),
                    reason="iflas-cupcake env not available (cdna_cupcake)")
def test_real_collapse_bridge_compatibility(tmp_path):
    """Compatibility coverage for genuine cupcake GFF/group bridge output.

    The main driver now reads pu_filter's hq bed directly; this keeps the helper
    checked against real cupcake ids for debugging/transition use.
    """
    from iflas.stages.align import minimap2_tgs
    from iflas.stages.collapse import collapse_isoforms

    d = str(tmp_path)
    syn.generate(d)
    bam = os.path.join(d, "flnc.bam")
    minimap2_tgs(f"{d}/genome.fa", f"{d}/flnc.fa", bam, platform="pacbio", threads=2)
    collapse_isoforms(f"{d}/flnc.fa", bam, f"{d}/tofu", fl_coverage=2)
    grp = os.path.join(d, "tofu.collapsed.group.txt")
    gff = os.path.join(d, "tofu.collapsed.good.gff")
    assert os.path.isfile(gff) and os.path.isfile(grp)

    bed = os.path.join(d, "iso.bed12+")
    collapsed_gff_to_bed12plus(gff, grp, bed)
    ref_gpe = os.path.join(d, "ref.gpe")
    _build_ref_gpe(os.path.join(d, "genes.gtf"), ref_gpe)

    out = find_as_lr(bed, grp, ref_gpe, d)
    assert "geneB:2200@2500-2650@3000" in [r[3] for r in _rows(out["SE"])]
    assert "geneD:300-700" in [r[3] for r in _rows(out["IR"])]
    assert "geneC:4000-4050" in [r[3] for r in _rows(out["A5SS"])]
    assert "geneE:4950-5000" in [r[3] for r in _rows(out["A3SS"])]
