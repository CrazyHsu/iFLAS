"""Collapse stage test on synthetic data: align -> cupcake collapse -> check isoform structures.

Requires the isolated `iflas-cupcake` conda env (cupcake routed via conda run).
"""

import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict

import pytest

from iflas.stages.align import minimap2_tgs
from iflas.stages.collapse import collapse_isoforms
from iflas.cli import main as cli_main

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import make_synthetic_data as syn  # noqa: E402


def _cupcake_available():
    if shutil.which("conda") is None:
        return False
    r = subprocess.run(["conda", "run", "-n", "iflas-cupcake", "collapse_isoforms_by_sam.py", "--help"],
                       capture_output=True)
    return r.returncode == 0


pytestmark = pytest.mark.skipif(not _cupcake_available(),
                                reason="iflas-cupcake env not available (cdna_cupcake)")


def _exons_per_transcript(gff):
    counts = defaultdict(int)
    with open(gff) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 9 and cols[2] == "exon":
                m = re.search(r'transcript_id "([^"]+)"', cols[8])
                if m:
                    counts[m.group(1)] += 1
    return counts


def test_collapse_keeps_distinct_isoforms(tmp_path):
    d = str(tmp_path)
    syn.generate(d)
    bam = os.path.join(d, "flnc.bam")
    minimap2_tgs(f"{d}/genome.fa", f"{d}/flnc.fa", bam, platform="pacbio", threads=2)
    gff = collapse_isoforms(f"{d}/flnc.fa", bam, f"{d}/tofu", fl_coverage=2)
    assert os.path.isfile(gff) and os.path.getsize(gff) > 0

    counts = _exons_per_transcript(gff)
    exon_counts = sorted(counts.values())
    # Designed isoforms collapse to distinct structures (>=7 core isoforms survive the dedup)
    assert len(counts) >= 7
    assert max(exon_counts) == 3   # geneA / geneB.iso1 are 3-exon
    assert 2 in exon_counts        # multiple 2-exon isoforms (A5SS/A3SS/geneD.iso1)
    assert min(exon_counts) == 1   # geneD.iso2_IR is single-exon (intron retained)


def test_cli_collapse_smoke_uses_configured_sample_layout(tmp_path):
    d = tmp_path / "data"
    d.mkdir()
    syn.generate(str(d))

    base = tmp_path / "out" / "demo" / "s1"
    mapping_dir = base / "mapping"
    mapping_dir.mkdir(parents=True)
    bam = mapping_dir / "flnc.mm2.sorted.bam"
    minimap2_tgs(str(d / "genome.fa"), str(d / "flnc.fa"), str(bam),
                 platform="pacbio", threads=2)

    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        f"""\
reference:
  - ref_strain: B73
    ref_genome: {d / "genome.fa"}
    ref_gtf: {d / "genes.gtf"}
collapse:
  fl_coverage: 2
dirs:
  out_dir: {tmp_path / "out"}
"""
    )
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample_name\tproject_name\tref_strain\tdata_processed_location\n"
        f"s1\tdemo\tB73\t{d / 'flnc.fa'}\n"
    )

    assert cli_main(["collapse", "-c", str(cfg), "-s", str(samples)]) == 0
    gff = base / "collapse" / "tofu.collapsed.good.gff"
    assert gff.is_file() and gff.stat().st_size > 0
