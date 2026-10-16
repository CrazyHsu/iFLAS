"""Acceptance tests for synthetic dataset v2 (scripts/make_synthetic_data_v2.py).

1. the generator is deterministic;
2. the truth tables are self-consistent with the files they describe;
3. the migrated core chain (init-ref -> mapping -> collapse -> refine -> pu_filter -> find_as) recovers
   exactly the designed AS events whose truth row says it should, and stays silent on the ones that
   document known limitations (R-2 mono-exon IR, offset < 5, MXE, < 2 reads, aligner microexons) or that
   the published mapping junction-combination filter removes (no:junc-filter, BUG-143);
4. the NGS confident filter drops alignment-wobble events; the PU filter never drops annotated isoforms;
5. the mapping junction filter keeps exactly the reads of the isoforms the truth says pass, and the NGS
   junction evidence it relied on agrees with the designed support on every "> 10 reads" decision.

Tests 3-4 need the full toolchain + the iflas-cupcake env (tool-gated like test_orchestrator_e2e.py).
"""
from __future__ import annotations

import csv
import hashlib
import os
import re
import shutil
import subprocess
import sys

import pandas as pd
import pytest
import yaml

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
V2 = os.path.join(ROOT, "test_data", "synthetic_v2")
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import make_synthetic_data_v2 as synv2  # noqa: E402

AS_FILES = {"SE": "SE.bed12+", "IR": "IR.bed6+", "A5SS": "A5SS.bed6+", "A3SS": "A3SS.bed6+"}


def _tsv(name):
    with open(os.path.join(V2, "truth", name), newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _rows(path):
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return [l.rstrip("\n").split("\t") for l in f if l.strip()]


# ----------------------------------------------------------------------------- 1. determinism
# at half depth some filter-deciding junctions fall into the generator's ambiguous NGS band (it warns; the committed
# full-depth data is asserted band-free in test_truth_tables_are_self_consistent)
@pytest.mark.filterwarnings("ignore:designed NGS support of filter-deciding junctions")
def test_generator_is_deterministic(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    synv2.generate(str(a), scale=0.5)
    synv2.generate(str(b), scale=0.5)
    for root, _, files in os.walk(a):
        for fn in files:
            pa = os.path.join(root, fn)
            pb = os.path.join(str(b), os.path.relpath(pa, str(a)))
            assert hashlib.sha256(open(pa, "rb").read()).hexdigest() == hashlib.sha256(open(pb, "rb").read()).hexdigest(), fn


# ----------------------------------------------------------------------------- 2. truth consistency
def test_truth_tables_are_self_consistent():
    isoforms = {r["isoform"]: r for r in _tsv("isoforms.tsv")}
    genes = {r["gene"] for r in _tsv("genes.tsv")}
    assert {r["gene"] for r in isoforms.values()} <= genes

    # genes.gtf holds exactly the annotated isoforms, with the designed exons
    gtf_exons = {}
    for row in _rows(os.path.join(V2, "genes.gtf")):
        name = row[8].split('transcript_id "')[1].split('"')[0]
        gtf_exons.setdefault(name, []).append((int(row[3]) - 1, int(row[4])))
    annotated = {n for n, r in isoforms.items() if r["annotated"] == "1"}
    assert set(gtf_exons) == annotated
    for name in annotated:
        assert ";".join(f"{s}-{e}" for s, e in sorted(gtf_exons[name])) == isoforms[name]["exons"]

    # every FLNC read belongs to a known isoform and the per-isoform counts match
    reads = _tsv("read_truth.tsv")
    counts = {}
    for r in reads:
        assert r["isoform"] in isoforms
        counts[r["isoform"]] = counts.get(r["isoform"], 0) + 1
    for name, r in isoforms.items():
        assert counts.get(name, 0) == int(r["n_flnc"]), name
    fasta_ids = {l[1:].strip() for l in open(os.path.join(V2, "flnc.fa")) if l.startswith(">")}
    assert fasta_ids == {r["read_id"] for r in reads}

    # events reference existing isoforms of the same gene and carry a decision for both outputs; an isoform
    # that the legacy mapping junction filter removes makes its events undetectable
    assert {r["junction_filter"] for r in isoforms.values()} <= {"pass", "fail:chain-support", "fail:consensus"}
    for ev in _tsv("as_events.tsv"):
        for iso in (ev["inc_isoforms"], ev["exc_isoforms"]):
            assert iso in isoforms and isoforms[iso]["gene"] == ev["gene"], ev["event_id"]
        assert ev["type"] in ("SE", "IR", "A5SS", "A3SS", "MXE")
        assert ev["expected_detection"] and ev["expected_confident"]
        assert int(ev["start"]) < int(ev["end"])
        filtered = any(isoforms[i]["junction_filter"] != "pass" for i in (ev["inc_isoforms"], ev["exc_isoforms"]))
        assert (ev["expected_detection"] == "no:junc-filter") == filtered, ev["event_id"]
    for row in _tsv("junction_support.tsv"):
        carriers = row["isoforms"].split(";")
        assert all(i in isoforms for i in carriers)
        assert row["consensus"] in ("0", "1")
        assert int(row["flnc_reads"]) == sum(int(isoforms[i]["n_flnc"]) for i in carriers)
        # a junction whose consensus hinges on NGS must sit well clear of the legacy "> 10 reads" decision
        if row["annotated"] == "0" and int(row["flnc_reads"]) < 3 and any(int(isoforms[i]["n_flnc"]) >= 2 for i in carriers):
            assert not (synv2.NGS_BAND[0] < int(row["ngs_reads_total"]) < synv2.NGS_BAND[1]), row
    for row in _tsv("diff_psi.tsv") + _tsv("allele_bias.tsv") + _tsv("palen.tsv"):
        assert row["expected_significant"] in ("0", "1")
    for row in _tsv("pa_sites.tsv"):
        assert 0 < float(row["frac"]) <= 1


# ----------------------------------------------------------------------------- 3/4. core chain
_TOOLS = all(shutil.which(t) for t in ("gtfToGenePred", "bedtools", "minimap2", "hisat2", "hisat2-build",
                                        "stringtie", "regtools", "samtools", "seqkit", "perl", "iflas"))


def _cupcake_available():
    if shutil.which("conda") is None:
        return False
    return subprocess.run(["conda", "run", "-n", "iflas-cupcake", "collapse_isoforms_by_sam.py", "--help"],
                          capture_output=True).returncode == 0


needs_chain = pytest.mark.skipif(not (_TOOLS and _cupcake_available()),
                                 reason="needs full core toolchain + iflas-cupcake env")


def _r_available():
    if shutil.which("conda") is None:
        return False
    return subprocess.run(["conda", "run", "-n", "iflas-diffsplice", "Rscript", "--version"],
                          capture_output=True).returncode == 0


needs_r = pytest.mark.skipif(not _r_available(), reason="needs the iflas-diffsplice env (Rscript)")


def _v2_config(work, filter_score, *, merge=True, data=V2):
    """`iflas init-ref` into work/ref and work/config.yaml from the committed synthetic config; returns (ref dir,
    sample-sheet rows with absolute paths). ``data`` is the synthetic v2 directory (or a copy of it)."""
    ref = work / "ref"
    ref.mkdir(parents=True, exist_ok=True)
    subprocess.run(["iflas", "init-ref", "-f", f"{data}/genome.fa", "--gtf", f"{data}/genes.gtf", "-o", str(ref), "--hisat2"],
                   check=True, capture_output=True)
    cfg = yaml.safe_load(open(os.path.join(V2, "config.yaml")))
    cfg["reference"][0].update({"ref_genome": f"{data}/genome.fa", "ref_gtf": f"{data}/genes.gtf",
                                "ref_size": f"{ref}/genome.len.txt", "ref_gpe": f"{ref}/genome.gpe", "ref_bed": f"{ref}/genome.bed",
                                "ref_mm2_index": f"{ref}/genome.mm2", "hisat2_index": f"{ref}/genome_hisat2"})
    cfg["dirs"]["out_dir"] = str(work / "iflas_out")
    cfg["pu_filter"]["filter_score"] = filter_score
    cfg["options"]["threads"] = 2
    cfg["options"]["merge_data_from_same_strain"] = merge
    with open(work / "config.yaml", "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    with open(os.path.join(V2, "samples.tsv"), newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    for row in rows:
        for k in ("data_processed_location", "primer", "ngs_left_reads", "ngs_right_reads"):
            if row[k] not in ("None", ""):
                row[k] = os.path.join(data, row[k])
    return ref, rows


def _write_samples(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def _iflas(*args):
    r = subprocess.run(["iflas", *args], capture_output=True, text=True)
    assert r.returncode == 0, (args[0], r.stderr[-4000:])
    return r


def _run_chain(work, filter_score):
    """init-ref + `iflas all` core stages on synthetic v2 with the committed config/sample sheet; returns the sample base dir."""
    ref, rows = _v2_config(work, filter_score)
    _write_samples(work / "samples.tsv", rows)
    _iflas("all", "-c", str(work / "config.yaml"), "-s", str(work / "samples.tsv"), "--stages", "mapping,collapse,refine,pu_filter,find_as", "-j", "1")
    return next((work / "iflas_out").glob("*/*/hqIsoforms")).parent, ref


def _rerun_pu_and_find_as(base, ref, dest, filter_score):
    """Re-run only pu_filter + find_as on a copy of a finished base dir (mapping/collapse/refine reused)."""
    shutil.copytree(base, dest)
    for sub in ("hqIsoforms", "as_events"):
        shutil.rmtree(dest / sub, ignore_errors=True)
    for cmd in (["iflas", "pu_filter", "-b", str(dest), "--ref-genome", f"{V2}/genome.fa", "--ref-bed", f"{ref}/genome.bed",
                 "--filter-score", str(filter_score)],
                ["iflas", "find_as", "-b", str(dest), "-r", f"{ref}/genome.gpe"]):
        r = subprocess.run(cmd, capture_output=True, text=True)
        assert r.returncode == 0, (cmd[1], r.stderr[-3000:])
    return dest


def _events(lr_dir, suffix=""):
    found = {}
    for typ, fn in AS_FILES.items():
        name = fn.replace(".bed", f".{suffix}bed") if suffix else fn
        for c in _rows(os.path.join(lr_dir, name)):
            found[(typ, c[0], int(c[1]), int(c[2]))] = c[3]
    return found


@pytest.fixture(scope="module")
def chain(tmp_path_factory):
    return _run_chain(tmp_path_factory.mktemp("v2_fs0"), 0.0)


@pytest.fixture(scope="module")
def chain_fs0(chain):
    return chain[0]


@pytest.fixture(scope="module")
def chain_fs05(chain, tmp_path_factory):
    base, ref = chain
    return _rerun_pu_and_find_as(base, ref, tmp_path_factory.mktemp("v2_fs05") / base.name, 0.5)


@needs_chain
def test_core_chain_matches_designed_events(chain_fs0):
    """With PU filtering disabled (filter_score=0) find_as must call exactly the designed events."""
    lr = chain_fs0 / "as_events" / "ordinary_as" / "LR"
    found, confident = _events(lr), _events(lr, "confident.")
    truth = _tsv("as_events.tsv")
    for ev in truth:
        key = (ev["type"], ev["chrom"], int(ev["start"]), int(ev["end"]))
        for got, col in ((found, "expected_detection"), (confident, "expected_confident")):
            want = ev[col] in ("yes", "both:A5SS+A3SS")
            assert (key in got) == want, f"{col}={ev[col]} but {'found' if key in got else 'missing'}: {ev['event_id']} ({ev['note']})"
    truth_keys = {(e["type"], e["chrom"], int(e["start"]), int(e["end"])) for e in truth}
    assert not (set(found) - truth_keys), f"events not designed: {set(found) - truth_keys}"


@needs_chain
def test_core_chain_hq_set_and_mono_exon_limitation(chain_fs0):
    """filter_score=0 keeps every multi-exon isoform; mono-exon isoforms never reach HQ (R-2)."""
    grouped = _rows(chain_fs0 / "refine" / "isoformGrouped.bed12+")
    hq = {r[3] for r in _rows(chain_fs0 / "hqIsoforms" / "hq.collapsed.bed12+")}
    multi = {r[3] for r in grouped if int(r[9]) > 1}
    mono = {r[3] for r in grouped if int(r[9]) == 1}
    assert multi <= hq | {h.split(",")[0] for h in hq}
    assert not (mono & hq)
    assert len(mono) == sum(1 for r in _tsv("isoforms.tsv") if r["exons"].count("-") == 1 and int(r["n_flnc"]) > 0)


@needs_chain
def test_pu_filter_keeps_annotated_and_reports_novel_recall(chain_fs05):
    """At the default threshold every annotated isoform survives; novel recall is reported (PU behaviour, BUG-107/141)."""
    feats = {r[0]: r for r in _rows(chain_fs05 / "hqIsoforms" / "isoFeatures.txt")[1:]}
    valid = {r[0] for r in _rows(chain_fs05 / "hqIsoforms" / "validIsoforms.lst")[1:]}
    annotated = {k for k, r in feats.items() if r[5] == "annotated"}
    novel = {k for k, r in feats.items() if r[5] == "novel"}
    assert annotated <= valid, sorted(annotated - valid)[:5]
    recall = len(novel & valid) / max(1, len(novel))
    print(f"\nPU novel recall at filter_score=0.5: {len(novel & valid)}/{len(novel)} = {recall:.2f}")
    assert 0.0 <= recall <= 1.0


@needs_chain
def test_mapping_junction_filter_matches_truth(chain_fs0):
    """mapping re-applies the published junction-combination filter (BUG-143): the kept FLNC reads are exactly
    the reads of the isoforms whose truth says they pass, and every designed junction's NGS evidence falls on the
    same side of the legacy `> 10 reads` threshold as the designed support."""
    kept = {l.strip() for l in open(chain_fs0 / "mapping" / "flnc.filterReads.lst") if l.strip()}
    status = {r["isoform"]: r["junction_filter"] for r in _tsv("isoforms.tsv")}
    expected = {r["read_id"] for r in _tsv("read_truth.tsv") if status[r["isoform"]] == "pass"}
    assert kept == expected, f"dropped: {sorted(expected - kept)[:5]} kept: {sorted(kept - expected)[:5]}"
    assert len(_rows(chain_fs0 / "mapping" / "flnc.addCVandID.bed12+")) == len(kept)
    assert (chain_fs0 / "mapping" / "flnc.mm2.sorted.bam.bai").is_file()

    assert len(list((chain_fs0 / "mapping" / "rna-seq" / "alignment").glob("repeat*.sorted.bam"))) == 1   # all 4 libraries in one run
    detected = {}
    for c in _rows(chain_fs0 / "mapping" / "rna-seq" / "reassembly" / "junctions.raw.bed"):
        sizes, starts = c[10].split(","), c[11].split(",")
        key = (c[0], int(c[1]) + int(sizes[0]), int(c[1]) + int(starts[1]))
        detected[key] = max(detected.get(key, 0), int(c[4]))       # the filter keeps a junction if any row > 10
    deltas, band = [], 0
    for r in _tsv("junction_support.tsv"):
        if r["ngs_detectable"] == "0":          # flanked by a microexon the aligner cannot place (R-11)
            continue
        key = (r["chrom"], int(r["start"]), int(r["end"]))
        designed, got = int(r["ngs_reads_total"]), detected.get(key, 0)
        deltas.append(designed - got)
        lo, hi = synv2.NGS_BAND
        if designed >= hi:
            assert got > synv2.NGS_JUNC_MIN_READS, (key, designed, got)
        elif designed <= lo:
            assert got <= synv2.NGS_JUNC_MIN_READS, (key, designed, got)
        else:
            band += 1                           # non-deciding junction near the threshold: no claim
    print(f"\ndesigned - detected NGS junction reads: min {min(deltas)} max {max(deltas)} "
          f"mean {sum(deltas) / len(deltas):.2f} over {len(deltas)} detectable junctions ({band} inside the band)")


@needs_chain
def test_refine_reads_level_assignment_matches_truth(chain_fs0):
    """refine's reads-level readsAssign (legacy "for reads" block, F-0): every read that survived the
    junction filter is assigned once. readsAssign gives a read its gene when its isoform is annotated, shares
    an intron with an annotated isoform of the gene, or is >= 90% covered by that gene's annotated exons;
    otherwise the read gets a novel cluster id (the known R-2/R-9 behaviour: unannotated mono-exon gene,
    mono-exon IR of a 2-exon gene, the second mutually-exclusive exon)."""
    kept = {l.strip() for l in open(chain_fs0 / "mapping" / "flnc.filterReads.lst") if l.strip()}
    isoforms = {r["isoform"]: r for r in _tsv("isoforms.tsv")}
    reads = {r["read_id"]: r["isoform"] for r in _tsv("read_truth.tsv")}
    exons = {n: [tuple(map(int, x.split("-"))) for x in r["exons"].split(";")] for n, r in isoforms.items()}
    annotated_introns, annotated_exons = {}, {}
    for name, r in isoforms.items():
        if r["annotated"] == "1":
            annotated_introns.setdefault(r["gene"], set()).update(r["junc_chain"].split(";"))
            annotated_exons.setdefault(r["gene"], []).extend(exons[name])

    def covered(iso):
        """fraction of the isoform's exonic bases inside the gene's annotated exons"""
        total = hit = 0
        for s, e in exons[iso["isoform"]]:
            total += e - s
            hit += sum(1 for p in range(s, e) if any(a <= p < b for a, b in annotated_exons.get(iso["gene"], [])))
        return hit / total

    rows = _rows(chain_fs0 / "refine" / "reads.assigned.unambi.bed12+")
    assert {row[3] for row in rows} == kept and len(rows) == len(kept)
    novel_isoforms = set()
    for row in rows:
        iso = isoforms[reads[row[3]]]
        shares = set(iso["junc_chain"].split(";")) & annotated_introns.get(iso["gene"], set())
        if iso["annotated"] == "1" or (iso["junc_chain"] and shares) or covered(iso) >= 0.9:
            assert row[14] == iso["gene"], (row[3], row[14])
        else:
            assert ":" in row[14], (row[3], row[14])        # novel gene cluster id
            novel_isoforms.add(iso["isoform"])
    assert novel_isoforms == {"monoC.full", "diA.ir", "mxeA.mxe2", "mxeB.mxe2"}
    gtf_tx = {l.split('transcript_id "')[1].split('"')[0]
              for l in open(chain_fs0 / "refine" / "reads.unambi.gtf") if "\ttranscript\t" in l}
    assert gtf_tx == kept


@pytest.fixture(scope="module")
def chain_pa(chain):
    """find_pa on the filter_score=0 chain (needs the R env)."""
    base, ref = chain
    if not _r_available():
        pytest.skip("needs the iflas-diffsplice env (Rscript)")
    r = subprocess.run(["iflas", "find_pa", "-b", str(base), "--ref-genome", f"{V2}/genome.fa",
                        "--ref-size", f"{ref}/genome.len.txt", "-t", "2"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]
    return base, ref


@needs_chain
@needs_r
def test_find_pa_recovers_designed_pa_sites(chain_pa):
    """`iflas find_pa -b` (legacy find_pa driver, F-3): every designed PA site with >= 8 supporting reads is
    called within 12 bp (TES jitter sd 3, PAC window 3), no PA site is called > 30 bp from a designed one,
    and the PA-group tables / motif profiles / plots are produced."""
    chain_fs0, _ = chain_pa
    pa_dir = chain_fs0 / "as_events" / "pa"
    called = [(c[0], c[5], int(c[2]), int(c[4])) for c in _rows(pa_dir / "PA.bed6+")]     # chrom, strand, site, reads
    designed = [(r["chrom"], r["strand"], int(r["site"]), int(r["n_reads_designed"]), r["gene"]) for r in _tsv("pa_sites.tsv")]
    for chrom, strand, site, n, gene in designed:
        if n >= 8:
            assert any(c == chrom and s == strand and abs(p - site) <= 12 for c, s, p, _ in called), (gene, site, n)
    for chrom, strand, site, n in called:
        assert any(c == chrom and s == strand and abs(p - site) <= 30 for c, s, p, _, _ in designed), (chrom, site, n)
    assert sum(1 for *_, n in called if n >= 5) == len(called)                              # PA.bed6+ applies the read filter
    for name in ("paCluster.bed8+", "isoform.paGrouped.tsv", "reads.paGrouped.tsv", "reads.paGrouped.bed6",
                 "reads.3endRevised.bed12+", "adjacentCleavageDistr.pdf"):
        assert (pa_dir / name).stat().st_size > 0, name
    assert len(_rows(pa_dir / "reads.3endRevised.bed12+")) == len(_rows(chain_fs0 / "refine" / "reads.assigned.unambi.bed12+"))
    for name in ("A.nucleotide", "AATAAA.PAS", "Other.PAS", "nucleotide.pdf", "PAS.pdf"):
        assert (pa_dir / "motif" / name).stat().st_size > 0, name
    apa_genes = {c[2] for c in _rows(pa_dir / "reads.paGrouped.tsv") if "," in c[3]}
    print(f"\nPA sites called: {len(called)}; designed with >= 8 reads: {sum(1 for d in designed if d[3] >= 8)}; APA genes: {sorted(apa_genes)}")


@needs_chain
@needs_r
def test_characterize_back_half_matches_truth(chain_pa):
    """`iflas characterize -b -r --ref-genome` (legacy charaterize_as back half, F-1): the AS events implied by
    the reference annotation are exactly the designed events whose two isoforms are both annotated; the
    confident LR events split into Known (reference) / Novel accordingly; all APA genes are Novel (every
    designed gene has one annotated TTS); splice-site motif tables, distance-to-TTS table and plots exist."""
    base, ref = chain_pa
    r = subprocess.run(["iflas", "characterize", "-b", str(base), "-r", f"{ref}/genome.gpe", "--ref-genome", f"{V2}/genome.fa"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]
    char = base / "as_events" / "characterization"
    isoforms = {r["isoform"]: r for r in _tsv("isoforms.tsv")}
    truth = _tsv("as_events.tsv")
    both_annotated = {(e["type"], e["chrom"], int(e["start"]), int(e["end"])) for e in truth
                      if e["type"] != "MXE" and all(isoforms[i]["annotated"] == "1" for i in (e["inc_isoforms"], e["exc_isoforms"]))}
    reference = {(t, c[0], int(c[1]), int(c[2])) for t in AS_FILES for c in _rows(char / f"{t}.reference.bed6+")}
    assert reference == both_annotated

    lr = base / "as_events" / "ordinary_as" / "LR"
    for t in AS_FILES:
        confident = {(t, c[0], int(c[1]), int(c[2])): c[3] for c in _rows(lr / AS_FILES[t].replace(".bed", ".confident.bed"))}
        if t == "IR":
            confident = {(t, c[0], int(c[1]), int(c[2])): c[3] for c in _rows(lr / "IR.bed6+")}      # IR has no NGS filter
        rows = [c for c in _rows(char / f"statistic.{t}.tsv") if not c[0].startswith("#")]
        assert len(rows) == len(confident), t
        known = {c[3] for c in rows if c[2] == "Known"}
        assert known == {name for key, name in confident.items() if key in both_annotated}, t
    apa = [c for c in _rows(char / "statistic.APA.tsv") if not c[0].startswith("#")]
    assert apa and all(c[2] == "Novel" for c in apa)
    assert {c[3] for c in apa} == {c[2] for c in _rows(base / "as_events" / "pa" / "reads.paGrouped.tsv") if "," in c[3]}
    assert len(_rows(char / "lrPA2TTS.tsv")) == len(_rows(base / "as_events" / "pa" / "reads.paGrouped.bed6"))
    for name in ("IR.tmp.txt", "SE.inc.tmp.txt", "SE.exc.tmp.txt", "A5SS.inc.tmp.txt", "A3SS.inc.tmp.txt"):
        fractions = [float(c[1]) for c in _rows(char / name)]
        assert fractions and abs(sum(fractions) - 1.0) < 1e-9, name
    assert _rows(char / "SE.inc.tmp.txt")[0][0] == "GT-AG"
    for pdf in ("SE.inc.ssMotif.pdf", "IR.ssMotif.pdf", "lrPA2TTS.pdf", "SE.pdf"):
        assert (char / pdf).stat().st_size > 0, pdf


@needs_chain
def test_visual_as_plots_target_genes_with_ngs_tracks(chain):
    """`iflas visual_as -b -r -g` (F-2, pyGenomeTracks): one PDF per target gene, region covering the designed gene
    span, the gene's HQ isoforms in the isoform track, NGS coverage + junction tracks present; `--merge` stacks
    the sample under the reference in isoViewer_sample_merged/."""
    base, ref = chain
    genes = ["seA", "irB", "a5B", "monoA", "antiA", "antiB", "nearA"]
    r = subprocess.run(["iflas", "visual_as", "-b", str(base), "-r", f"{ref}/genome.gpe", "-g", ",".join(genes), "-t", "2"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]
    span = {g["gene"]: (g["chrom"], int(g["start"]), int(g["end"])) for g in _tsv("genes.tsv")}
    hq_by_gene = {}
    for c in _rows(base / "hqIsoforms" / "hq.collapsed.bed12+"):
        hq_by_gene.setdefault(c[13], []).append(c[3])          # col 14 = gene name (13 = cupcake gene id)
    for gene in genes:
        d = base / "visual_as" / gene
        assert (d / f"{gene}.pdf").stat().st_size > 0, gene
        ini = (d / "tracks.ini").read_text()
        assert "[reference]" in ini and "[coverage_repeat0_" in ini, gene
        assert ("[junctions_" in ini) == (gene != "monoA"), gene          # no NGS junction inside a mono-exon gene
        assert {c[3] for c in _rows(d / f"{gene}.iflas_isoforms.bed")} == set(hq_by_gene.get(gene, [])), gene   # monoA: no HQ isoform (R-2)
        chrom, start, end = span[gene]
        assert ini.count("file_type = bedgraph") == 1                       # the merged sample maps one NGS repeat
        ref_rows = _rows(d / f"{gene}.reference.bed")
        assert ref_rows[0][0] == chrom and min(int(x[1]) for x in ref_rows) <= start and max(int(x[2]) for x in ref_rows) >= end

    r = subprocess.run(["iflas", "visual_as", "--merge", "-b", str(base), "-r", f"{ref}/genome.gpe", "-g", "seA,irB"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]
    merged = base.parent.parent / "isoViewer_sample_merged"
    for gene in ("seA", "irB"):
        assert (merged / gene / f"{gene}.pdf").stat().st_size > 0
        assert f"[isoforms_{base.name}]" in (merged / gene / "tracks.ini").read_text()


@pytest.fixture(scope="module")
def chain_conditions(tmp_path_factory):
    """Per-sample run (no strain merge): mapping for all four samples (three are NGS-only), the long-read
    chain for the TGS sample condA_rep1, then `iflas diff_as` on condA vs condB."""
    if not _r_available():
        pytest.skip("needs the iflas-diffsplice env (rMATS)")
    work = tmp_path_factory.mktemp("v2_cond")
    ref, rows = _v2_config(work, 0.0, merge=False)
    _write_samples(work / "samples.tsv", rows)
    _write_samples(work / "samples_tgs.tsv", [r for r in rows if r["tgs_plat"] != "None"])
    for samples, stages in (("samples.tsv", "mapping"), ("samples_tgs.tsv", "collapse,refine,pu_filter,find_as")):
        _iflas("all", "-c", str(work / "config.yaml"), "-s", str(work / samples), "--stages", stages, "-j", "2")
    # --pu-filter: merge the HQ isoforms. The default (isoforms of *confident* AS events) is filtered per sample
    # with that sample's own NGS, and condA_rep1's library is inclusion-dominant by design, so the
    # condition-specific skip junctions (<= 10 reads in condA) are not NGS-confident there (R-13).
    _iflas("diff_as", "-c", str(work / "config.yaml"), "-s", str(work / "samples.tsv"), "-d", "condA,condB", "--pu-filter", "-t", "2")
    return work


@needs_chain
def test_diff_as_calls_designed_differential_events(chain_conditions):
    """F-4: rMATS on condA vs condB (2 x 2 NGS replicates) over the merged HQ isoform GTF of the TGS sample. Every
    designed strong event (|dPSI| >= 0.4) that rMATS can define from the GTF (fromGTF.<type>.txt) is significant,
    no designed null gene is, and at least two strong events are called. rMATS cannot define a two-exon skip (seC) as
    SE, nor an A5SS/A3SS whose shared flanking exon is a terminal exon with isoform-specific ends (a5B; R-14)."""
    das = chain_conditions / "iflas_out" / "das"
    comp = das / "synv2_condA_vs_synv2_condB"
    assert (das / "all_sample_merged_iso.gtf").stat().st_size > 0 and (comp / "SE.MATS.JC.txt").is_file()
    rmats_type = {"IR": "RI", "SE": "SE", "A3SS": "A3SS", "A5SS": "A5SS"}
    defined = {t: {str(g).strip('"') for g in pd.read_csv(comp / f"fromGTF.{rmats_type[t]}.txt", sep="\t")["GeneID"]} for t in rmats_type}
    sig = {t: {str(g).strip('"') for g in pd.read_csv(das / "synv2_condA_vs_synv2_condB.sigDiffAS" / f"{t}.sig.txt", sep="\t")["GeneID"]}
           for t in rmats_type}
    truth = _tsv("diff_psi.tsv")
    strong = {e["gene"]: e["event_id"].split(":")[1] for e in truth if e["expected_significant"] == "1" and abs(float(e["delta"])) >= 0.4}
    null = {e["gene"] for e in truth if e["expected_significant"] == "0"}
    called = set().union(*sig.values())
    undefined = []
    for gene, typ in strong.items():
        types = (typ,) if typ in ("SE", "IR") else ("A5SS", "A3SS")
        if not any(gene in defined[t] for t in types):
            undefined.append((gene, typ))
            continue
        assert any(gene in sig[t] for t in types), (gene, typ, {t: sorted(v) for t, v in sig.items()})
    assert len(strong) - len(undefined) >= 2, (strong, undefined)
    assert not (null & called), null & called
    weak = {e["gene"] for e in truth if e["expected_significant"] == "1" and abs(float(e["delta"])) < 0.4}
    print(f"\nrMATS significant genes: { {t: sorted(v) for t, v in sig.items()} }; strong events rMATS could not define: {undefined}; "
          f"weak designed (|dPSI|<0.4) called: {sorted(weak & called)}/{sorted(weak)}")


def _cupcake_phaser_available():
    if shutil.which("conda") is None:
        return False
    return subprocess.run(["conda", "run", "-n", "iflas-cupcake", "run_phaser.py", "--help"], capture_output=True).returncode == 0


@needs_chain
@pytest.mark.skipif(not _cupcake_phaser_available(), reason="needs IsoPhase (run_phaser.py) in the iflas-cupcake env")
def test_asas_phases_designed_allele_biased_events(chain):
    """F-5: `iflas asas -b --ref-genome --flnc` phases the FL reads of every locus with >= 10 reads (IsoPhase) and
    relates haplotypes to the SE isoform pairs: the four designed allele-biased genes (hap bias 0.8-0.85) are
    reported with their two haplotypes, the two designed null genes (bias 0.5) are not."""
    base, ref = chain
    r = subprocess.run(["iflas", "asas", "-b", str(base), "--ref-genome", f"{V2}/genome.fa", "--flnc", f"{V2}/flnc.fa", "-t", "2"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]
    loci = sorted(p.name for p in (base / "asas" / "by_loci").glob("*size*"))
    assert len(loci) >= 6, loci
    phased = [p for p in (base / "asas" / "by_loci").glob("*size*") if (p / "phased.partial.cleaned.human_readable.txt").is_file()]
    rows = _rows(base / "asas" / "partialAsRelatedHaplotype.txt")
    reported = {r[2] for r in rows}                      # refGene column
    expected = {r["gene"] for r in _tsv("allele_bias.tsv") if r["expected_significant"] == "1"}
    null = {r["gene"] for r in _tsv("allele_bias.tsv") if r["expected_significant"] == "0"}
    print(f"\nIsoPhase loci: {len(loci)} (phased with SNPs: {len(phased)}); ASAS genes reported: {sorted(reported)}")
    assert not (null & reported), null & reported
    assert expected <= reported, (expected - reported, rows)
    for row in rows:
        assert row[1] == "SE" and row[4] != row[6] and row[5] and row[7], row

    # F-5b: allele-specific expression from the merged sample's NGS reads (one repeat holding all four libraries)
    with open(os.path.join(V2, "samples.tsv"), newline="") as f:
        libs = [r for r in csv.DictReader(f, delimiter="\t")]
    left = ",".join(os.path.join(V2, r["ngs_left_reads"]) for r in libs)
    right = ",".join(os.path.join(V2, r["ngs_right_reads"]) for r in libs)
    r = subprocess.run(["iflas", "asas", "-b", str(base), "--ase", "--ref-genome", f"{V2}/genome.fa", "--ngs-left", left, "--ngs-right", right, "-t", "2"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]
    merged = _rows(base / "ase" / "repeat0.mergedASEfile.txt")
    by_gene = {row[1]: row for row in merged}
    ase_truth = {r["gene"]: r for r in _tsv("ase.tsv")}
    assert set(by_gene) == set(ase_truth), sorted(by_gene)
    for gene, row in by_gene.items():
        tgs_ref, tgs_alt, ngs_ref, ngs_alt = float(row[2]), float(row[3]), float(row[5]), float(row[7])
        assert tgs_ref + tgs_alt > 0 and ngs_ref + ngs_alt > 0, row
        tgs_share, ngs_share = max(tgs_ref, tgs_alt) / (tgs_ref + tgs_alt), max(ngs_ref, ngs_alt) / (ngs_ref + ngs_alt)
        imbalanced = ase_truth[gene]["expected_imbalanced"] == "1"
        # the ASAS genes put each isoform on its own haplotype, so their total expression per haplotype is balanced;
        # the ASE gene (both isoforms 85% hapA) is imbalanced in both the FL reads and the salmon NGS counts
        assert (tgs_share >= 0.7) == imbalanced, (gene, row)
        assert (ngs_share >= 0.7) == imbalanced, (gene, row)       # salmon counts on ~60 reads wobble a few percent
        if imbalanced:
            assert (tgs_ref > tgs_alt) == (ngs_ref > ngs_alt), row                # same haplotype dominates in both
    print(f"ASE table: { {g: r[2:] for g, r in sorted(by_gene.items())} }")


@needs_chain
def test_palen_as_calls_designed_tail_length_shifts(chain_pa):
    """F-6: `iflas palen_as -b --palen polya_results.tsv -pa_sup 5` tests the poly(A) tail lengths between the
    isoforms of every confident AS event with the published rule (t-test and Kruskal both <= 0.001 on PASS reads,
    >= pa_sup reads each side). The designed shifts that the rule can resolve at the designed depth
    (truth palen.tsv `detectable_published`) are called, nothing else is."""
    base, _ = chain_pa
    r = subprocess.run(["iflas", "palen_as", "-b", str(base), "--palen", f"{V2}/polya_results.tsv", "-pa_sup", "5"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]
    sig_dir = base / "palenAS" / "mergeByJunc"
    called = {}
    for t in AS_FILES:
        for c in _rows(sig_dir / f"{t}.palenAndAS.sig.bed12+"):
            called.setdefault(t, set()).add(c[13])              # gene name of the representative HQ isoform (col 14)
    truth = _tsv("palen.tsv")
    expected = {(e["event_id"].split(":")[1], e["gene"]) for e in truth if e["detectable_published"] == "1"}
    assert {(t, g) for t, genes in called.items() for g in genes} == expected, (called, expected)
    designed = {e["gene"] for e in truth if e["expected_significant"] == "1"}
    print(f"\npalen_as: called {sorted(expected)}; designed shifts below the published rule's reach: {sorted(designed - {g for _, g in expected})}")
    assert (base / "palenAS" / "palenAPA" / "apaRelatedPalen.txt").is_file()


def _pacbio_tools_available():
    if shutil.which("conda") is None:
        return False
    return subprocess.run(["conda", "run", "-n", "iflas-pacbio", "lima", "--version"], capture_output=True).returncode == 0


@needs_chain
@pytest.mark.skipif(not _pacbio_tools_available(), reason="needs the iflas-pacbio env (lima / isoseq)")
def test_preproc_pacbio_recovers_designed_flnc_from_ccs_bam(tmp_path):
    """F-7: `iflas preproc` on the synthetic CCS BAM (already-called reads: ccs skipped, lima --isoseq, isoseq refine
    --require-polya) yields exactly the designed FLNC set -- reads carrying both primers and a >=20 nt poly(A).
    Sequences equal flnc.fa except that isoseq cannot tell genomic 3' A's from the tail and trims them too
    (BUG-101/102 closed against the real lima/isoseq)."""
    cfg = yaml.safe_load(open(os.path.join(V2, "config.yaml")))
    cfg["reference"][0].update({"ref_genome": f"{V2}/genome.fa", "ref_gtf": f"{V2}/genes.gtf"})
    cfg["dirs"]["out_dir"] = str(tmp_path / "iflas_out")
    with open(tmp_path / "config.yaml", "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    with open(os.path.join(V2, "samples.tsv"), newline="") as f:
        row = next(r for r in csv.DictReader(f, delimiter="\t") if r["tgs_plat"] == "pacbio")
    ngs = tmp_path / "ngs"                              # a legacy fastp copy beside the inputs would be reused: isolate
    shutil.copytree(os.path.join(V2, "ngs"), ngs)
    row.update({"data_location": f"{V2}/ccs.bam", "data_processed_location": "None", "primer": f"{V2}/primer.fa",
                "ngs_left_reads": str(ngs / os.path.basename(row["ngs_left_reads"])),
                "ngs_right_reads": str(ngs / os.path.basename(row["ngs_right_reads"])), "use_fmlrc2": "false"})
    _write_samples(tmp_path / "samples.tsv", [row])
    r = subprocess.run(["iflas", "preproc", "-c", str(tmp_path / "config.yaml"), "-s", str(tmp_path / "samples.tsv"), "-t", "2"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]
    pacbio = tmp_path / "iflas_out" / "synv2" / row["sample_name"] / "preprocess" / "pacbio"
    assert (pacbio / "fl.5p--3p.bam").is_file() and (pacbio / "flnc.bam").is_file()
    flnc = {}
    with open(pacbio / "rawFlnc.fq") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                name = line[1:].split()[0]
            elif i % 4 == 1:
                flnc[name] = line.strip()
    designed = {f"m64000_000000_000000/{t['zmw']}/ccs": t for t in _tsv("read_truth.tsv") if t["flnc_recovered"] == "1"}
    assert set(flnc) == set(designed), (len(flnc), len(designed), sorted(set(flnc) ^ set(designed))[:3])
    seqs = {}
    with open(f"{V2}/flnc.fa") as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].strip()
            else:
                seqs[name] = line.strip()
    identical = trimmed = 0
    bad = []
    for n, t in designed.items():
        got, want = flnc[n], seqs[t["read_id"]]
        if got == want:
            identical += 1
        elif want.startswith(got) and set(want[len(got):]) <= {"A"}:
            trimmed += 1                                   # genomic A's adjacent to the tail went with the poly(A)
        else:
            bad.append(n)
    assert not bad, (len(bad), bad[:3])
    assert trimmed == sum(1 for t in designed.values() if seqs[t["read_id"]].endswith("A"))
    print(f"\npreproc PacBio: {len(flnc)} FLNC reads == designed (both primers, poly(A) >= 20); "
          f"{identical} sequences identical to flnc.fa, {trimmed} shorter only by genomic 3' A's")


# ----------------------------------------------------------------------------- F-8. the whole chain in one command
@pytest.fixture(scope="module")
def chain_all(tmp_path_factory):
    """Default `iflas all` (no --stages: preproc ... palen_as) on the committed merged config, run where paths are
    hostile: the inputs are a copy of synthetic v2 and the outputs a work dir, both under names with a space and
    Chinese characters (every shell command quotes its paths). The copied NGS reads also keep fastp from reusing a
    product beside the repo data; the TGS row gets the synthetic poly(A) table so palen_as has input."""
    if not (_r_available() and _cupcake_phaser_available()):
        pytest.skip("needs the iflas-diffsplice (R) and iflas-cupcake (IsoPhase) envs")
    root = tmp_path_factory.mktemp("v2_all")
    data, work = root / "input data 输入", root / "run dir 运行"
    shutil.copytree(V2, data, ignore=shutil.ignore_patterns("truth"))
    work.mkdir()
    ref, rows = _v2_config(work, 0.0, data=data)
    for row in rows:
        row["polya_location"] = f"{data}/polya_results.tsv" if row["tgs_plat"] != "None" else "None"
    _write_samples(work / "samples.tsv", rows)
    r = _iflas("all", "-c", str(work / "config.yaml"), "-s", str(work / "samples.tsv"), "-j", "1")
    return work, r.stdout


@needs_chain
def test_iflas_all_runs_the_whole_published_chain(chain_all, chain_fs0):
    """F-8: one default `iflas all` runs preproc per sample row, then every migrated stage on the merged sample.
    diff_as is skipped with a notice (merged mapping pools the four libraries, so no per-condition alignments);
    the AS events equal the stage-by-stage core chain's."""
    work, stdout = chain_all
    base = work / "iflas_out" / "synv2" / "SYN_all"
    assert len(list((work / "iflas_out" / "synv2").glob("*/preprocess/rna-seq/fastp.*_R1.fastq.gz"))) == 4   # preproc ran per row (BUG-152 layout)
    for rel in ("mapping/flnc.mm2.sorted.bam", "collapse/tofu.collapsed.good.gff", "refine/isoformGrouped.bed12+",
                "refine/reads.assigned.unambi.bed12+", "hqIsoforms/hq.collapsed.bed12+", "as_events/ordinary_as/LR/SE.bed12+",
                "as_events/pa/PA.bed6+", "as_events/characterization/lrPA2TTS.tsv", "asas/partialAsRelatedHaplotype.txt",
                "ase/repeat0.mergedASEfile.txt"):
        assert (base / rel).is_file() and (base / rel).stat().st_size > 0, rel
    assert (base / "palenAS" / "palenAPA" / "apaRelatedPalen.txt").is_file()      # no APA tail-length shift is designed
    assert any((base / "palenAS" / "mergeByJunc").glob("*.palenAndAS.sig.bed12+"))
    assert any((base / "visual_as").glob("*/*.pdf"))
    lines = {l.split(":")[0].strip(): l for l in stdout.splitlines() if l.startswith("  ")}
    assert "skipped: per-sample NGS alignments not found" in lines["diff_as"], lines["diff_as"]
    for stage in ("asas", "ase", "palen_as", "find_pa", "characterize", "visual_as"):
        assert "skipped" not in lines[stage], lines[stage]
    assert _events(base / "as_events" / "ordinary_as" / "LR") == _events(chain_fs0 / "as_events" / "ordinary_as" / "LR")
    print("\niflas all: preproc + 11 per-sample stages ran, diff_as skipped with a notice (merged mapping); AS events == core chain")


def _numeric_chromosome_copy(dest):
    """synthetic v2 with Ensembl-style chromosome names ("chr1" -> "1", as Ensembl Plants maize, rice and Arabidopsis
    references name them); only the genome and the annotation name chromosomes."""
    shutil.copytree(V2, dest, ignore=shutil.ignore_patterns("truth"))
    for name in ("genome.fa", "genome.unmasked.fa", "genes.gtf"):
        path = dest / name
        path.write_text(re.sub(r"^(>?)chr(\d+)", r"\1\2", path.read_text(), flags=re.M))
    return dest


@pytest.fixture(scope="module")
def chain_numeric(tmp_path_factory):
    """Default `iflas all` and `iflas report -all` on the numerically named copy (otherwise as chain_all)."""
    if not (_r_available() and _cupcake_phaser_available()):
        pytest.skip("needs the iflas-diffsplice (R) and iflas-cupcake (IsoPhase) envs")
    root = tmp_path_factory.mktemp("v2_numeric")
    data, work = _numeric_chromosome_copy(root / "data"), root / "run"
    work.mkdir()
    ref, rows = _v2_config(work, 0.0, data=data)
    for row in rows:
        row["polya_location"] = f"{data}/polya_results.tsv" if row["tgs_plat"] != "None" else "None"
    _write_samples(work / "samples.tsv", rows)
    _iflas("all", "-c", str(work / "config.yaml"), "-s", str(work / "samples.tsv"), "-j", "1")
    _iflas("report", "-c", str(work / "config.yaml"), "-s", str(work / "samples.tsv"), "-all")
    return work / "iflas_out"


@needs_chain
def test_numeric_chromosome_names_give_the_same_results(chain_all, chain_numeric):
    """All synthetic data named chromosomes "chr1"...; with Ensembl-style "1".."4" the whole chain and the report plots
    run, and the products equal the "chr" run's once the names are mapped back."""
    work, _ = chain_all
    chr_base, num_base = work / "iflas_out" / "synv2" / "SYN_all", chain_numeric / "synv2" / "SYN_all"
    rels = ["hqIsoforms/hq.collapsed.bed12+", "refine/isoformGrouped.bed12+", "as_events/pa/PA.bed6+",
            "as_events/characterization/lrPA2TTS.tsv", "asas/partialAsRelatedHaplotype.txt", "ase/repeat0.mergedASEfile.txt"]
    rels += [f"as_events/ordinary_as/LR/{fn}" for fn in AS_FILES.values()]
    rels += [str(p.relative_to(chr_base)) for p in (chr_base / "palenAS" / "mergeByJunc").glob("*.palenAndAS.sig.bed12+")]
    for rel in rels:
        renamed = sorted(re.sub(r"(?<![A-Za-z0-9])chr(\d+)(?!\d)", r"\1", line) for line in open(chr_base / rel))
        assert renamed == sorted(open(num_base / rel)), rel
    sample = next(p for p in (chain_numeric / "reports").iterdir() if (p / "basicStatistics").is_dir())
    for folder, n in {"basicStatistics": 5, "asPattern": 2, "geneStrucPlots": 1, "alleleAsPlots": 1, "palenAsPlots": 1}.items():
        assert len([p for p in (sample / folder).glob("*.pdf") if p.stat().st_size > 0]) >= n, folder


@needs_chain
def test_pu_filter_model_selection_on_synthetic_v2(chain, tmp_path):
    """F-10: synthetic v2 has too few reliable negatives for the legacy model selection (every designed novel isoform
    carries >= 5 % of its gene's reads and its junctions are NGS-supported), so `--select-best-model` falls back to
    GB with a warning instead of crashing in train_test_split; on a feature table with ten weak novel isoforms the
    real selection runs: five classifiers, a chosen model, per-model thresholds and the ROC plot."""
    base, ref = chain
    dest = tmp_path / base.name
    shutil.copytree(base, dest)
    shutil.rmtree(dest / "hqIsoforms")
    r = _iflas("pu_filter", "-b", str(dest), "--ref-genome", f"{V2}/genome.fa", "--ref-bed", f"{ref}/genome.bed",
               "--select-best-model", "--auto-filter-score", "--draw-auc")
    assert "WARNING: select_best_model needs" in r.stdout and "model GB" in r.stdout
    assert len(_rows(dest / "hqIsoforms" / "hq.collapsed.bed12+")) > 0

    from iflas.pu.iso_pu import selectBestModel
    feats = pd.read_csv(dest / "hqIsoforms" / "isoFeatures.txt", sep="\t", index_col=0)
    weak = feats.index[feats.annotation == "novel"][:10]
    feats.loc[weak, "ratioIsoToGene"] = 0.01                     # ten reliable negatives
    work = tmp_path / "select"
    work.mkdir()
    cwd = os.getcwd()
    os.chdir(work)
    try:
        model, thresholds = selectBestModel(feats, drawAUC=True, selectPUscore=True)
    finally:
        os.chdir(cwd)
    assert model in ("RF", "GB", "DT", "SVM", "NB") and set(thresholds) == {"RF", "GB", "DT", "SVM", "NB"}
    assert all(len(v) == 5 and all(0 <= t <= 1 for t in v) for v in thresholds.values())
    assert (work / "PU.bagging.ROC_AUC.pdf").stat().st_size > 0
    print(f"\nselect_best_model on synthetic v2 features (10 weak negatives): {model}; thresholds "
          + ", ".join(f"{m}={sum(v) / len(v):.2f}" for m, v in thresholds.items()))


@needs_chain
@pytest.mark.skipif(not _pacbio_tools_available(), reason="needs the iflas-pacbio env (ccs / lima / isoseq)")
def test_preproc_runs_ccs_on_fabricated_subreads(tmp_path):
    """U2: `iflas preproc` on a fabricated Sequel subreads BAM (scripts/make_fake_pacbio_subreads.py: one insert per
    ZMW read 6-10 times with 3 % errors, PacBio tags incl. kinetics) runs pbindex -> dataset create -> ccs -> lima ->
    isoseq refine; the consensus FLNC sequences must reproduce the designed inserts (poly(A)-trimmed, allowing isoseq
    to trim genomic 3' A's)."""
    import make_fake_pacbio_subreads as fake
    subreads = tmp_path / "fake.subreads.bam"
    fake.main(["--flnc", f"{V2}/flnc.fa", "--out", str(subreads), "--zmws", "30"])
    cfg = yaml.safe_load(open(os.path.join(V2, "config.yaml")))
    cfg["reference"][0].update({"ref_genome": f"{V2}/genome.fa", "ref_gtf": f"{V2}/genes.gtf"})
    cfg["dirs"]["out_dir"] = str(tmp_path / "iflas_out")
    cfg["ccs"] = {"min_pass": 2, "min_ccs_length": 50, "min_predicted_accuracy": 0.9}
    with open(tmp_path / "config.yaml", "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    with open(os.path.join(V2, "samples.tsv"), newline="") as f:
        row = next(r for r in csv.DictReader(f, delimiter="\t") if r["tgs_plat"] == "pacbio")
    row.update({"data_location": str(subreads), "data_processed_location": "None", "primer": f"{V2}/primer.fa",
                "ngs_left_reads": "None", "ngs_right_reads": "None", "use_fmlrc2": "false"})
    _write_samples(tmp_path / "samples.tsv", [row])
    r = _iflas("preproc", "-c", str(tmp_path / "config.yaml"), "-s", str(tmp_path / "samples.tsv"), "-t", "2")
    pacbio = tmp_path / "iflas_out" / "synv2" / row["sample_name"] / "preprocess" / "pacbio"
    assert (pacbio / "CCS.bam").is_file() and (pacbio / "fl.5p--3p.bam").is_file() and (pacbio / "flnc.bam").is_file()
    flnc = {}
    with open(pacbio / "rawFlnc.fq") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                name = line[1:].split()[0]
            elif i % 4 == 1:
                flnc[name.split("/")[1]] = line.strip()                           # keyed by ZMW
    truth = {r["zmw"]: r["flnc_read"] for r in csv.DictReader(open(str(subreads) + ".truth.tsv"), delimiter="\t")}
    seqs, name = {}, None
    for line in open(f"{V2}/flnc.fa"):
        if line.startswith(">"):
            name = line[1:].strip()
        else:
            seqs[name] = line.strip()
    assert len(flnc) == len(truth) == 30, (len(flnc), len(truth))
    exact = sum(1 for zmw, read in truth.items() if seqs[read].startswith(flnc[zmw]) and set(seqs[read][len(flnc[zmw]):]) <= {"A"})
    assert exact >= 27, exact                                                     # ccs consensus of 6-10 noisy passes
    print(f"\nccs on fabricated subreads: {len(flnc)} FLNC reads, {exact}/30 identical to the designed insert (modulo 3' A's)")


@needs_chain
@pytest.mark.skipif(shutil.which("nanopolish") is None, reason="needs nanopolish (core iflas env)")
def test_palen_as_calls_nanopolish_on_simulated_squiggles(tmp_path):
    """U2: `polya_len_calling` (nanopolish index + polya through the ToolRunner) on a fabricated multi-read fast5 whose
    signal is simulated from nanopolish's own RNA pore model (scripts/make_fake_fast5.py: leader, adapter, poly(A)
    of a designed length, transcript 3'->5'): the designed tail lengths are recovered."""
    import make_fake_fast5 as fake
    from scipy.stats import spearmanr
    from iflas.stages.palen_as import polya_len_calling, read_polya_lengths
    ont = tmp_path / "ont"
    fake.main(["--flnc", f"{V2}/flnc.fa", "--model", f"{V2}/r9.4_70bps.u_to_t_rna.5mer.template.model", "--out-dir", str(ont), "--reads", "20"])
    bam = ont / "reads.bam"
    subprocess.run(f"minimap2 -ax splice -uf -k14 {V2}/genome.fa {ont}/reads.fastq 2>/dev/null | samtools sort -o {bam} - && samtools index {bam}",
                   shell=True, check=True)
    out = polya_len_calling(str(ont / "fast5"), str(ont / "reads.fastq"), str(bam), f"{V2}/genome.fa", tmp_path / "palen",
                            sequencing_summary=str(ont / "sequencing_summary.txt"), threads=2)
    rows = _rows(out)
    assert rows[0][:2] == ["readname", "contig"] and "qc_tag" in rows[0] and len(rows) == 21
    tags = rows[0]
    passed = {r[0]: float(r[tags.index("polya_length")]) for r in rows[1:] if r[tags.index("qc_tag")] == "PASS"}
    truth = {r["read_id"]: int(r["polya_len"]) for r in csv.DictReader(open(ont / "tail_truth.tsv"), delimiter="\t")}
    assert len(passed) >= 10, len(passed)
    est, des = zip(*[(passed[k], truth[k]) for k in passed])
    rho, mae = spearmanr(est, des).correlation, sum(abs(a - b) for a, b in zip(est, des)) / len(est)
    assert rho >= 0.9 and mae <= 6, (rho, mae)
    assert read_polya_lengths([out]) == passed                                   # the downstream reader keeps PASS rows
    print(f"\nnanopolish on simulated squiggles: {len(passed)}/20 PASS, Spearman {rho:.3f}, mean |error| {mae:.1f} nt")


def _report_config(dest, out_dir, ref):
    """config.yaml + samples.tsv in ``dest`` that point `iflas report` at ``out_dir``; returns the ``-c/-s`` arguments."""
    cfg = yaml.safe_load(open(os.path.join(V2, "config.yaml")))
    cfg["reference"][0].update({"ref_genome": f"{V2}/genome.fa", "ref_gtf": f"{V2}/genes.gtf", "ref_gpe": f"{ref}/genome.gpe", "ref_bed": f"{ref}/genome.bed"})
    cfg["dirs"]["out_dir"] = str(out_dir)
    with open(dest / "config.yaml", "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    with open(os.path.join(V2, "samples.tsv"), newline="") as f:
        rows = [r for r in csv.DictReader(f, delimiter="\t")]
    for row in rows:
        for k in ("data_processed_location", "primer", "ngs_left_reads", "ngs_right_reads"):
            if row[k] not in ("None", ""):
                row[k] = os.path.join(V2, row[k])
    _write_samples(dest / "samples.tsv", rows)
    return "-c", str(dest / "config.yaml"), "-s", str(dest / "samples.tsv")


@needs_chain
@needs_r
def test_report_novel_hq_as_plots_the_isoform_score_written_by_pu_filter(chain, tmp_path):
    """F-9: pu_filter now writes isoIdentity/isoformScore.txt, so `iflas report --novel-hq-as` produces its plot."""
    base, ref = chain
    assert len(_rows(base / "isoIdentity" / "isoformScore.txt")) == len(_rows(base / "hqIsoforms" / "isoFeatures.txt")) - 1
    r = _iflas("report", *_report_config(tmp_path, base.parent.parent, ref), "--novel-hq-as")
    assert "Generated report plots" in r.stdout
    pdfs = list((base.parent.parent / "reports").glob("*/isoformScore.pdf"))
    assert pdfs and pdfs[0].stat().st_size > 0


@pytest.fixture(scope="module")
def chain_reports(chain_pa, tmp_path_factory):
    """`iflas report -all` with real R on a copy of the synthetic chain. The unit tests replace R, which hid
    BUG-164..166; here the stages a report reads (characterize, visual_as, palen_as, asas) run on the copy first. Synthetic v2
    designs no tail-length difference between the PA sites of a gene, so one APA line (seB's two designed sites)
    is written to drive the APA plot."""
    if not _cupcake_phaser_available():
        pytest.skip("needs IsoPhase (run_phaser.py) in the iflas-cupcake env")
    base, ref = chain_pa
    work = tmp_path_factory.mktemp("v2_reports")
    shutil.copytree(base.parent.parent, work / "iflas_out", symlinks=True)
    copy = work / "iflas_out" / base.parent.name / base.name
    _iflas("characterize", "-b", str(copy), "-r", f"{ref}/genome.gpe", "--ref-genome", f"{V2}/genome.fa")
    _iflas("visual_as", "-b", str(copy), "-r", f"{ref}/genome.gpe", "-g", "seA,irB", "-t", "2")
    _iflas("palen_as", "-b", str(copy), "--palen", f"{V2}/polya_results.tsv", "-pa_sup", "5")
    _iflas("asas", "-b", str(copy), "--ref-genome", f"{V2}/genome.fa", "--flnc", f"{V2}/flnc.fa", "-t", "2")
    tails = ",".join(str(60 + i) for i in range(6)) + ";" + ",".join(str(110 + i) for i in range(6))
    (copy / "palenAS" / "palenAPA" / "apaRelatedPalen.txt").write_text(f"seB\t-\tchr1_6050;chr1_5980\t6;6\t{tails}\n")
    r = _iflas("report", *_report_config(work, work / "iflas_out", ref), "-all")
    return work / "iflas_out" / "reports", r.stdout


@needs_chain
@needs_r
def test_report_all_renders_every_plot_with_real_r(chain_reports):
    """Every report helper writes non-empty PDFs with the real R stack (ggplot2 4, txdbmaker), and the HTML
    report collects them."""
    reports, stdout = chain_reports
    sample = next(p for p in reports.iterdir() if (p / "basicStatistics").is_dir())
    expected = {"basicStatistics": 5, "asPattern": 2, "geneStrucPlots": 1, "alleleAsPlots": 1, "palenAsPlots": 1, "palenApaPlots": 1}
    for folder, n in expected.items():
        pdfs = [p for p in (sample / folder).glob("*.pdf") if p.stat().st_size > 0]
        assert len(pdfs) >= n, (folder, sorted(p.name for p in (sample / folder).glob("*")))
    assert (sample / "isoformScore.pdf").stat().st_size > 0
    # the HTML pages show the basic statistics / AS pattern plots and one structure image per visual_as gene
    pngs = {p.name for p in (reports / "html").rglob("*.png") if p.stat().st_size > 0}
    wanted = {f"SYN_all_basic_{p.stem}.png" for folder in ("basicStatistics", "asPattern") for p in (sample / folder).glob("*.pdf")}
    wanted |= {f"SYN_all_gene_{g}_isoformStruc.png" for g in ("seA", "irB")}
    assert wanted <= pngs, wanted - pngs
    assert (reports / "html" / "iflas_report.html").is_file() and "Wrote HTML report" in stdout


@needs_chain
@needs_r
def test_report_diff_renders_with_real_r(chain_conditions, tmp_path):
    """`iflas report -diff` plots the significant differential AS events of condA vs condB with real R."""
    _iflas("report", "-c", str(chain_conditions / "config.yaml"), "-s", str(chain_conditions / "samples.tsv"), "-diff")
    pdfs = list((chain_conditions / "iflas_out" / "das").glob("*.sigDiffAS/sigDiff.AS_distribution.pdf"))
    assert pdfs and pdfs[0].stat().st_size > 0


@pytest.fixture(scope="module")
def chain_single_end(tmp_path_factory):
    """Default `iflas all` on the TGS sample alone, its NGS given as two single-end repeats (legacy `;` syntax):
    hisat2 -U per repeat, `samtools cat` + `stringtie --merge` over two repeats, salmon -r per repeat."""
    if not (_r_available() and _cupcake_phaser_available()):
        pytest.skip("needs the iflas-diffsplice (R) and iflas-cupcake (IsoPhase) envs")
    work = tmp_path_factory.mktemp("v2_single_end")
    ref, rows = _v2_config(work, 0.0, merge=False)
    row = next(r for r in rows if r["tgs_plat"] != "None")
    row["ngs_left_reads"] = f"{V2}/ngs/condA_rep1_R1.fastq.gz;{V2}/ngs/condA_rep2_R1.fastq.gz"
    row["ngs_right_reads"], row["ngs_reads_paired"], row["polya_location"] = "None", "single", f"{V2}/polya_results.tsv"
    _write_samples(work / "samples.tsv", [row])
    r = _iflas("all", "-c", str(work / "config.yaml"), "-s", str(work / "samples.tsv"), "-j", "1")
    return work / "iflas_out" / "synv2" / row["sample_name"], r.stdout


@needs_chain
def test_single_end_ngs_in_two_repeats_runs_the_whole_chain(chain_single_end):
    """Every stage runs on single-end NGS split in two repeats, and allele-tagged NGS reads are found again
    (BUG-173: HISAT2 keeps the /1 of single-end names, so the ASE read lists matched nothing and every count was 0)."""
    base, stdout = chain_single_end
    lines = {l.split(":")[0].strip(): l for l in stdout.splitlines() if l.startswith("  ")}
    assert "skipped" in lines["diff_as"] and all("skipped" not in lines[s] for s in ("mapping", "pu_filter", "find_as", "ase", "palen_as"))
    rna = base / "mapping" / "rna-seq"
    for rel in ("alignment/repeat0.sorted.bam", "alignment/repeat1.sorted.bam", "alignment/repeat0.gtf", "alignment/repeat1.gtf",
                "reassembly/stringtie_merged.gtf", "reassembly/junctions.bed"):
        assert (rna / rel).stat().st_size > 0, rel
    assert _rows(base / "as_events" / "ordinary_as" / "NGS" / "SE.bed12+")                 # NGS-confirmed events exist
    ase_truth = {r["gene"]: r for r in _tsv("ase.tsv")}
    for i in (0, 1):
        by_gene = {row[1]: row for row in _rows(base / "ase" / f"repeat{i}.mergedASEfile.txt")}
        assert all(float(row[5]) + float(row[7]) > 0 for row in by_gene.values()), by_gene
        for gene, t in ase_truth.items():
            ngs_ref, ngs_alt = float(by_gene[gene][5]), float(by_gene[gene][7])
            assert (max(ngs_ref, ngs_alt) / (ngs_ref + ngs_alt) >= 0.7) == (t["expected_imbalanced"] == "1"), (i, by_gene[gene])


def _synthetic_gene2go(path):
    """gene -> GO table over the synthetic genes: generic terms for all, and RNA splicing (GO:0008380) for the designed
    differential-splicing genes plus six others (enricher's default minGSSize is 10)."""
    genes = sorted({r["gene"] for r in _tsv("genes.tsv")})
    diff = sorted({r["gene"] for r in _tsv("diff_psi.tsv") if r["expected_significant"] == "1"})
    splicing = set(diff) | set([g for g in genes if g not in diff][:6])
    rows = [(g, t) for i, g in enumerate(genes) for t in ("GO:0005737", "GO:0008150") + (("GO:0005634",) if i % 2 else ())]
    rows += [(g, "GO:0008380") for g in sorted(splicing)]
    path.write_text("gene\tgo\n" + "".join(f"{g}\t{t}\n" for g, t in rows))
    return path


@needs_chain
@needs_r
def test_go_enrichment_with_real_r(chain_conditions, tmp_path):
    """`iflas go` (one and two gene lists) and `iflas diff_as --go` run clusterProfiler for real: RNA splicing is
    enriched among the significant differential-splicing genes (BUG-170..172 were only visible this way)."""
    gene2go = _synthetic_gene2go(tmp_path / "gene2go.tsv")
    (tmp_path / "a.lst").write_text("a3A\nirB\nseA\n")
    (tmp_path / "b.lst").write_text("seA\nseC\ncpxA\na5B\n")
    _iflas("go", "-tg", str(tmp_path / "a.lst"), "-s", "one", "-bg", str(gene2go), "-o", "one", "--outdir", str(tmp_path))
    _iflas("go", "-tg", f"{tmp_path / 'a.lst'},{tmp_path / 'b.lst'}", "-s", "A,B", "-bg", str(gene2go), "-o", "two", "--outdir", str(tmp_path))
    one = pd.read_csv(tmp_path / "GO" / "one.goEnrichResults.txt", sep="\t")
    two = pd.read_csv(tmp_path / "GO" / "two.goEnrichResults.txt", sep="\t")
    assert one.loc[one.ID == "GO:0008380", "p.adjust"].item() <= 0.05
    assert set(two.loc[two.ID == "GO:0008380", "Cluster"]) == {"A", "B"}                   # both lists kept (BUG-172)
    for name in ("one", "two"):
        assert (tmp_path / "GO" / f"{name}.goEnrichResults.pdf").stat().st_size > 0      # not ".pdf108" (BUG-170)
    _iflas("diff_as", "-c", str(chain_conditions / "config.yaml"), "-s", str(chain_conditions / "samples.tsv"), "-d", "condA,condB",
           "--pu-filter", "-t", "2", "--go", "--gene2go", str(gene2go))
    sig = next((chain_conditions / "iflas_out" / "das").glob("*.sigDiffAS"))
    das_go = pd.read_csv(sig / "GO" / "sigDiff.goEnrichResults.txt", sep="\t")
    assert "GO:0008380" in set(das_go.ID) and (sig / "GO" / "sigDiff.goEnrichResults.pdf").stat().st_size > 0


# ----------------------------------------------------------------------------- ONT direct RNA (synthetic_v2/ont)
ONT = os.path.join(V2, "ont")


def _ont_row(template, name, data_location, *, fmlrc=False, ngs=True):
    row = dict(template)
    row.update({"sample_name": name, "replicate": name, "tgs_plat": "nanopore", "strategy": "None", "primer": "None",
                "data_processed_location": "None", "data_location": data_location, "use_fmlrc2": "true" if fmlrc else "false",
                "polya_location": "None"})
    if not ngs:
        row["ngs_left_reads"] = row["ngs_right_reads"] = "None"
    return row


def _set_option(work, key, value):
    cfg = yaml.safe_load(open(work / "config.yaml"))
    cfg["options"][key] = value
    with open(work / "config.yaml", "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)


@pytest.fixture(scope="module")
def chain_ont_mixed(tmp_path_factory):
    """Default `iflas all` on one strain sequenced twice: the PacBio FLNC row (with NGS and the synthetic poly(A) table)
    and an ONT RNA004 row given as dorado's unaligned BAM (no NGS, tails in pt:i), poly(A) trimmed by default (BUG-178)."""
    if not (_r_available() and _cupcake_phaser_available()):
        pytest.skip("needs the iflas-diffsplice (R) and iflas-cupcake (IsoPhase) envs")
    work = tmp_path_factory.mktemp("v2_ont_mixed")
    ref, rows = _v2_config(work, 0.0, merge=True)
    pac = dict(next(r for r in rows if r["tgs_plat"] != "None"), polya_location=f"{V2}/polya_results.tsv")
    _write_samples(work / "samples.tsv", [pac, _ont_row(pac, "ont004", f"{ONT}/rna004/dorado.bam", ngs=False)])
    r = _iflas("all", "-c", str(work / "config.yaml"), "-s", str(work / "samples.tsv"), "-j", "1")
    return work / "iflas_out" / "synv2", r.stdout


@needs_chain
def test_ont_rna004_dorado_bam_merges_with_pacbio(chain_ont_mixed):
    """F-11: a dorado BAM is read (samtools fastq) and its pt:i tails become palen_as input; the strain merges an
    ONT row without NGS (BUG-176) and FASTQ with FASTA reads (BUG-177); trimmed tails keep ONT reads through collapse."""
    project, stdout = chain_ont_mixed
    ont = project / "ont004" / "preprocess" / "nanopore"
    assert open(ont / "rawFlnc.fq").read(1) == "@" and (ont / "flnc.polyA_trimmed.fq").stat().st_size > 0
    truth = {r["read_id"]: int(r["dorado_pt"]) for r in csv.DictReader(open(f"{ONT}/rna004/read_truth.tsv"), delimiter="\t")}
    table = pd.read_csv(ont / "dorado_polya.tsv", sep="\t")
    assert set(table.loc[table.qc_tag == "PASS", "readname"]) == {k for k, v in truth.items() if v > 0}
    base = project / "SYN_all"
    grouped = {r for line in open(base / "collapse" / "tofu.collapsed.group.txt") for r in line.rstrip("\n").split("\t")[1].split(",")}
    assert len(grouped & set(truth)) >= 0.95 * len(truth) and any("/ccs/" in r for r in grouped)     # both platforms collapsed
    found = _events(base / "as_events" / "ordinary_as" / "LR")
    truth_ev = {(e["type"], e["chrom"], int(e["start"]), int(e["end"])): e for e in _tsv("as_events.tsv")}
    want = {k for k, e in truth_ev.items() if e["expected_detection"] in ("yes", "both:A5SS+A3SS")}
    assert len(want & set(found)) >= len(want) - 2 and len(set(found) - set(truth_ev)) <= 2, (len(want & set(found)), set(found) - set(truth_ev))
    lines = {l.split(":")[0].strip(): l for l in stdout.splitlines() if l.startswith("  ")}
    assert "skipped" not in lines["palen_as"] and any((base / "palenAS" / "mergeByJunc").glob("*.palenAndAS.sig.bed12+"))


@pytest.fixture(scope="module")
def chain_ont_rna002(tmp_path_factory):
    """Default `iflas all` on an ONT RNA002 sample (guppy-style FASTQ with U) corrected by fmlrc2 with its NGS, with
    simulated fast5 signal for the reads of the genes whose tail lengths are designed to differ (nanopolish)."""
    if not (_r_available() and _cupcake_phaser_available() and shutil.which("nanopolish") and shutil.which("fmlrc2")):
        pytest.skip("needs R, IsoPhase, nanopolish and fmlrc2")
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    import make_fake_fast5
    work = tmp_path_factory.mktemp("v2_ont_rna002")
    ref, rows = _v2_config(work, 0.0, merge=False)
    pac = next(r for r in rows if r["tgs_plat"] != "None")
    genes = {r["gene"] for r in _tsv("palen.tsv") if r["expected_significant"] == "1"}
    seqs = {}
    with open(f"{ONT}/rna002/reads.fastq") as f:
        for i, line in enumerate(f):
            name = line[1:].strip() if i % 4 == 0 else name
            if i % 4 == 1:
                seqs[name] = line.strip()
    truth = [r for r in csv.DictReader(open(f"{ONT}/rna002/read_truth.tsv"), delimiter="\t") if r["gene"] in genes]
    records = [(r["read_id"], seqs[r["read_id"]][:len(seqs[r["read_id"]]) - int(r["emitted_tail"])], int(r["polya_len"])) for r in truth]
    row = _ont_row(pac, "ont002", f"{ONT}/rna002/reads.fastq", fmlrc=True)
    row["fast5_dir"] = make_fake_fast5.write_fast5(records, f"{V2}/r9.4_70bps.u_to_t_rna.5mer.template.model", str(work / "signal"))
    _write_samples(work / "samples.tsv", [row])
    r = _iflas("all", "-c", str(work / "config.yaml"), "-s", str(work / "samples.tsv"), "-j", "1")
    return work / "iflas_out" / "synv2" / "ont002", {rec[0] for rec in records}, r.stdout


@needs_chain
def test_ont_rna002_fmlrc2_and_nanopolish_tails(chain_ont_rna002):
    """fmlrc2 corrects the ONT reads with the NGS (ropebwt2 -> fmlrc2-convert -> fmlrc2), and nanopolish reads the
    simulated signal against the raw reads' own alignments (BUG-175: the corrected reads' BAM cost tails as
    SUFFCLIP); the designed tail-length shifts the published rule can resolve are called."""
    base, signal_reads, stdout = chain_ont_rna002
    assert (base / "preprocess" / "fmlrc" / "fmlrc_corrected.fasta").stat().st_size > 0
    assert (base / "mapping" / "rawFlnc.mm2.sorted.bam").is_file()
    polya = pd.read_csv(base / "palenAS" / "polya_results.tsv", sep="\t")
    passed = set(polya.loc[polya.qc_tag == "PASS", "readname"])
    assert passed <= signal_reads and len(passed) >= 0.85 * len(signal_reads) and "SUFFCLIP" not in set(polya.qc_tag)
    called = set()
    for t in AS_FILES:
        called |= {(t, c[13]) for c in _rows(base / "palenAS" / "mergeByJunc" / f"{t}.palenAndAS.sig.bed12+")}
    assert {(e["event_id"].split(":")[1], e["gene"]) for e in _tsv("palen.tsv") if e["detectable_published"] == "1"} <= called


def test_preproc_trims_the_polya_of_gzipped_ont_reads(tmp_path):
    """Nanopore FASTQ given gzipped is decompressed (BUG-169) and, by default (BUG-178), a tail-trimmed copy
    becomes the processed reads: the cut lands within a few nt of the tail the basecaller emitted."""
    ref, rows = _v2_config(tmp_path, 0.0, merge=False)
    pac = next(r for r in rows if r["tgs_plat"] != "None")
    _write_samples(tmp_path / "samples.tsv", [_ont_row(pac, "ont004", f"{ONT}/rna004/reads.fastq.gz", ngs=False)])
    _iflas("all", "-c", str(tmp_path / "config.yaml"), "-s", str(tmp_path / "samples.tsv"), "--stages", "preproc")
    ont = tmp_path / "iflas_out" / "synv2" / "ont004" / "preprocess" / "nanopore"
    raw = open(ont / "rawFlnc.fq").read().splitlines()
    cut = open(ont / "flnc.polyA_trimmed.fq").read().splitlines()
    emitted = {r["read_id"]: int(r["emitted_tail"]) for r in csv.DictReader(open(f"{ONT}/rna004/read_truth.tsv"), delimiter="\t")}
    assert len(raw) == len(cut) == 4 * len(emitted)
    over = sorted(len(raw[i + 1]) - len(cut[i + 1]) - emitted[raw[i][1:]] for i in range(0, len(raw), 4))
    assert over[0] >= 0 and over[len(over) // 2] <= 5            # never short of the tail; median a few genomic A's more


def test_preproc_keeps_ont_tails_when_trimming_is_off(tmp_path):
    """options.ont_trim_polya: false is the published pipeline: the basecalled reads themselves are the processed reads."""
    ref, rows = _v2_config(tmp_path, 0.0, merge=False)
    pac = next(r for r in rows if r["tgs_plat"] != "None")
    _write_samples(tmp_path / "samples.tsv", [_ont_row(pac, "ont004", f"{ONT}/rna004/reads.fastq.gz", ngs=False)])
    _set_option(tmp_path, "ont_trim_polya", False)
    _iflas("all", "-c", str(tmp_path / "config.yaml"), "-s", str(tmp_path / "samples.tsv"), "--stages", "preproc")
    ont = tmp_path / "iflas_out" / "synv2" / "ont004" / "preprocess" / "nanopore"
    assert (ont / "rawFlnc.fq").stat().st_size > 0 and not (ont / "flnc.polyA_trimmed.fq").exists()
