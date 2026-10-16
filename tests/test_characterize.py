"""Unit tests for the migrated charaterize_as pure-core file helpers.

Covers the file-I/O characterization core ported line-by-line from legacy
``charaterize_as.py``:

  * ``iflas.io.files.filterFile``  (originally commonFuncs.py:76)
  * ``iflas.stages.characterize.getAnnoASList``
  * ``iflas.stages.characterize.getASstatistics``
  * ``iflas.stages.characterize.getSpliceSite``

External-tool branches (SE/A5SS/A3SS spliceSite decompose Perl helpers) are
checked via a recording fake runner for command contract only; the heavy
``charaterize_as`` driver (pybedtools + R + Perl + full ordinary_as/pa inputs)
remains a deferred real-data step.
"""

import os
import subprocess

import pytest

from iflas.io.files import filterFile
from iflas.stages.characterize import (
    build_anno_as_lists,
    characterize_sample,
    drawSSmotif,
    getAnnoASList,
    getASstatistics,
    getSpliceSite,
    splicesite2seq,
)
import shutil


def _read(path):
    with open(path) as f:
        return [line.rstrip("\n") for line in f]


# ---------------------------------------------------------------------------
# filterFile (io/files.py)
# ---------------------------------------------------------------------------

def test_filterfile_file_against_dict_intersection(tmp_path):
    origin = tmp_path / "origin.lst"
    origin.write_text("a\nb\nc\n")
    target = {"b": "row-b", "c": "row-c", "d": "row-d"}
    returned = filterFile(originFile=str(origin), targetFile=target,
                          outFile=str(tmp_path / "out.txt"), returnFlag=True)
    assert sorted(returned) == ["b", "c"]


def test_filterfile_mode_exclude(tmp_path):
    origin = tmp_path / "origin.lst"
    origin.write_text("a\nb\n")
    target = ["a", "b", "c", "d"]
    out = tmp_path / "excluded.txt"
    returned = filterFile(originFile=str(origin), targetFile=target, mode="e",
                          outFile=str(out), returnFlag=True)
    assert sorted(returned) == ["c", "d"]
    assert sorted(_read(str(out))) == ["c", "d"]


def test_filterfile_file_against_file_writes_target_rows(tmp_path):
    origin = tmp_path / "origin.lst"
    origin.write_text("k1\nk2\n")
    target = tmp_path / "target.tsv"
    target.write_text("k1\tval1\nk2\tval2\nk3\tval3\n")
    out = tmp_path / "out.txt"
    filterFile(originFile=str(origin), targetFile=str(target), outFile=str(out))
    assert sorted(_read(str(out))) == ["k1\tval1", "k2\tval2"]


def test_filterfile_writes_entries_in_target_order(tmp_path):
    """BUG-162: legacy wrote set order, which Python 3 changes from run to run; entries now follow the target."""
    origin = tmp_path / "origin.lst"
    origin.write_text("".join("ev{}\n".format(i) for i in range(0, 40, 2)))
    target = ["ev{}".format(i) for i in range(39, -1, -1)] + ["ev0"]
    kept = filterFile(originFile=str(origin), targetFile=target, outFile=str(tmp_path / "i.txt"), returnFlag=True)
    dropped = filterFile(originFile=str(origin), targetFile=target, mode="e", outFile=str(tmp_path / "e.txt"), returnFlag=True)
    assert kept == _read(str(tmp_path / "i.txt")) == ["ev{}".format(i) for i in range(38, -1, -2)]
    assert dropped == _read(str(tmp_path / "e.txt")) == ["ev{}".format(i) for i in range(39, 0, -2)]


# ---------------------------------------------------------------------------
# getAnnoASList
# ---------------------------------------------------------------------------

def test_get_anno_as_list_non_pa_skips_novel_and_builds_event(tmp_path):
    inFile = tmp_path / "IR.bed6+"
    inFile.write_text(
        "chr1\t100\t200\tgeneA:100-150\t.\t+\tgeneA\n"
        "chr2\t300\t400\tgeneB:300-350\t.\t-\tgeneB\tNovel\n"  # contains "Novel" -> skip
    )
    out = tmp_path / "LR.IR.lst"
    getAnnoASList(str(inFile), str(out))
    # asEvent = col0:col3:col5
    assert _read(str(out)) == ["chr1:geneA:100-150:+"]


def test_get_anno_as_list_uniq_dedup(tmp_path):
    inFile = tmp_path / "A5SS.bed6+"
    inFile.write_text(
        "chr1\t100\t200\tgeneA:100-150\t.\t+\tgeneA\n"
        "chr1\t100\t200\tgeneA:100-150\t.\t+\tgeneA\n"
    )
    out = tmp_path / "out.lst"
    getAnnoASList(str(inFile), str(out), uniq=True)
    assert _read(str(out)) == ["chr1:geneA:100-150:+"]


def test_get_anno_as_list_pa_requires_multiple_sites(tmp_path):
    inFile = tmp_path / "reads.paGrouped.tsv"
    inFile.write_text(
        "chr3\t+\tgeneC\t100,200\n"   # >1 PA site -> kept
        "chr3\t+\tgeneD\t500\n"       # single PA site -> skipped
    )
    out = tmp_path / "LR.APA.lst"
    getAnnoASList(str(inFile), str(out), PA=True)
    # asEvent = ":".join(cols[0:4])
    assert _read(str(out)) == ["chr3:+:geneC:100,200"]


def test_get_anno_as_list_append_preserves_existing_rows(tmp_path):
    """BUG-117: append=True should append LR events after known events."""
    out = tmp_path / "NGS.A5SS.lst"
    first = tmp_path / "known.bed6+"
    first.write_text("chr1\t10\t20\tgeneA:10-15\t.\t+\tgeneA\n")
    getAnnoASList(str(first), str(out))
    assert _read(str(out)) == ["chr1:geneA:10-15:+"]

    second = tmp_path / "lr.bed6+"
    second.write_text("chr2\t30\t40\tgeneB:30-35\t.\t-\tgeneB\n")
    getAnnoASList(str(second), str(out), append=True, uniq=True)
    assert _read(str(out)) == ["chr1:geneA:10-15:+", "chr2:geneB:30-35:-"]


# ---------------------------------------------------------------------------
# getASstatistics
# ---------------------------------------------------------------------------

def test_get_as_statistics_ir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    asFile = tmp_path / "IR.bed6+"
    asFile.write_text(
        "chr1\t100\t200\tgeneA:100-150\t.\t+\tgeneA\n"
        "chr2\t300\t400\tgeneB:300-350\t.\t-\tgeneB\n"
    )
    annoFile = tmp_path / "IR.anno.lst"
    annoFile.write_text("chr1:geneA:100-150:+\n")
    novelFile = tmp_path / "IR.novel.lst"
    novelFile.write_text("chr2:geneB:300-350:-\n")
    out = tmp_path / "statistic.IR.tsv"
    getASstatistics(asType="IR", asFile=str(asFile), annoFile=str(annoFile),
                    novelFile=str(novelFile), outFile=str(out))
    rows = _read(str(out))
    assert rows[0].startswith("##AS ID is composed of")
    assert "#Chr\tStrand\tKnown or Novel\tAS ID\tGene" in rows
    assert "chr1\t+\tKnown\tgeneA:100-150\tgeneA" in rows
    assert "chr2\t-\tNovel\tgeneB:300-350\tgeneB" in rows
    # the scratch filtered.txt must be cleaned up
    assert not (tmp_path / "filtered.txt").exists()


def test_get_as_statistics_se_uses_gene_col_13(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    asFile = tmp_path / "SE.confident.bed12+"
    asFile.write_text(
        "chr2\t1000\t2000\tgeneB:1100@1200-1300@1400\t.\t+\t"
        "1000\t2000\t0,0,0\t2\t100,100\t0,900\tgeneB\n"
    )
    annoFile = tmp_path / "SE.anno.lst"
    annoFile.write_text("chr2:geneB:1100@1200-1300@1400:+\n")
    novelFile = tmp_path / "SE.novel.lst"
    novelFile.write_text("")
    out = tmp_path / "statistic.SE.tsv"
    getASstatistics(asType="SE", asFile=str(asFile), annoFile=str(annoFile),
                    novelFile=str(novelFile), outFile=str(out))
    rows = _read(str(out))
    assert "chr2\t+\tKnown\tgeneB:1100@1200-1300@1400\tgeneB" in rows


def test_get_as_statistics_pa(tmp_path):
    annoFile = tmp_path / "paGroup.anno.tsv"
    annoFile.write_text(
        "chr3\t+\tgeneC\t100,200,300\n"
        "chr3\t+\tgeneD\t500\n"  # single PA -> skip
    )
    novelFile = tmp_path / "paGroup.novel.tsv"
    novelFile.write_text("chr4\t-\tgeneE\t700,800\n")
    out = tmp_path / "statistic.APA.tsv"
    getASstatistics(asType="PA", annoFile=str(annoFile), novelFile=str(novelFile),
                    outFile=str(out))
    rows = _read(str(out))
    assert rows[0] == "#Chr\tStrand\tKnown or Novel\tGene\tPA Sites"
    assert "chr3\t+\tKnown\tgeneC\t100,200,300" in rows
    assert "chr4\t-\tNovel\tgeneE\t700,800" in rows
    assert all("geneD" not in r for r in rows)


# ---------------------------------------------------------------------------
# getSpliceSite
# ---------------------------------------------------------------------------

def test_get_splice_site_ir_is_pure_transform(tmp_path):
    asFile = tmp_path / "LR.IR.lst"
    asFile.write_text("chr1:geneA:100-150:+\nchr2:geneB:300-350:-\n")
    out = tmp_path / "IR.splicesite"
    getSpliceSite(asType="IR", asFile=str(asFile), outFile=str(out))
    assert _read(str(out)) == ["chr1\t100\t150\t+", "chr2\t300\t350\t-"]


def test_get_splice_site_se_command_contract(tmp_path, monkeypatch):
    monkeypatch.setenv("IFLAS_UTILS_DIR", str(tmp_path))
    monkeypatch.chdir(tmp_path)
    (tmp_path / "confident.SE.lst").write_text("dummy\n")

    class RecordingRunner:
        def __init__(self):
            self.calls = []

        def run(self, cmd, **kwargs):
            self.calls.append((cmd, kwargs))
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    runner = RecordingRunner()
    getSpliceSite(asType="SE", asFile="confident.SE.lst", runner=runner)
    assert len(runner.calls) == 1
    cmd, _ = runner.calls[0]
    assert os.path.basename(cmd[-2] if cmd[-1] == "confident.SE.lst" else cmd[1]) == "seDecompose.pl"
    assert cmd[-1] == "confident.SE.lst"


def test_get_splice_site_a3ss_a5ss_use_matching_files_and_flags(tmp_path, monkeypatch):
    monkeypatch.setenv("IFLAS_UTILS_DIR", str(tmp_path))
    monkeypatch.chdir(tmp_path)
    (tmp_path / "custom.A3SS.lst").write_text("chr1:geneA:100-120:+\n")
    (tmp_path / "custom.A5SS.lst").write_text("chr1:geneA:200-220:+\n")

    class RecordingRunner:
        def __init__(self):
            self.calls = []

        def run(self, cmd, **kwargs):
            self.calls.append((cmd, kwargs))
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    runner = RecordingRunner()
    getSpliceSite(asType="A3SS", asFile="custom.A3SS.lst", runner=runner)
    getSpliceSite(asType="A5SS", asFile="custom.A5SS.lst", runner=runner)

    a3_cmd, a3_kwargs = runner.calls[0]
    assert os.path.basename(a3_cmd[1]) == "anssDecompose.pl"
    assert a3_cmd[2:4] == ["-n", "3"]
    assert a3_cmd[-1] == "custom.A3SS.lst"
    assert os.path.basename(a3_kwargs["stdout"].name) == "A3SS.inc.splicesite"
    assert os.path.basename(a3_kwargs["stderr"].name) == "A3SS.exc.splicesite"

    a5_cmd, a5_kwargs = runner.calls[1]
    assert os.path.basename(a5_cmd[1]) == "anssDecompose.pl"
    assert a5_cmd[2:4] == ["-n", "5"]
    assert a5_cmd[-1] == "custom.A5SS.lst"
    assert os.path.basename(a5_kwargs["stdout"].name) == "A5SS.inc.splicesite"
    assert os.path.basename(a5_kwargs["stderr"].name) == "A5SS.exc.splicesite"


# ---------------------------------------------------------------------------
# build_anno_as_lists  (front-half driver of charaterize_as)
# ---------------------------------------------------------------------------

class _RecordingRunner:
    def __init__(self):
        self.calls = []

    def run(self, cmd, **kwargs):
        self.calls.append((cmd if isinstance(cmd, str) else list(cmd), kwargs))
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")


def _write_ordinary_as_lr(base):
    """Write a minimal ordinary_as/LR tree (find_as output contract)."""
    lrDir = base / "ordinary_as" / "LR"
    lrDir.mkdir(parents=True)
    (lrDir / "SE.bed12+").write_text(
        "chr1\t1000\t2000\tgeneA:se1\t.\t+\textra\n"
        "chr1\t1000\t2000\tgeneA:seNovel\t.\t+\tNovel\n"  # contains "Novel" -> skip
    )
    (lrDir / "A5SS.bed6+").write_text(
        "chr1\t100\t200\tgeneA:a5shared\t.\t+\tgeneA\n"
        "chr1\t300\t400\tgeneA:a5lronly\t.\t+\tgeneA\n"
    )
    (lrDir / "A3SS.bed6+").write_text(
        "chr2\t500\t600\tgeneB:a3one\t.\t-\tgeneB\n"
    )
    (lrDir / "IR.bed6+").write_text(
        "chr3\t10\t20\tgeneC:ir1\t.\t+\tgeneC\n"
    )
    return base / "ordinary_as"


def test_build_anno_as_lists_no_ngs(tmp_path):
    ordinary = _write_ordinary_as_lr(tmp_path)
    out = tmp_path / "characterization"
    prev = os.getcwd()
    build_anno_as_lists(str(ordinary), str(out), ngs=False)

    # working directory restored after the driver finishes
    assert os.getcwd() == prev

    assert _read(str(out / "LR.SE.lst")) == ["chr1:geneA:se1:+"]
    assert _read(str(out / "LR.A5SS.lst")) == [
        "chr1:geneA:a5shared:+", "chr1:geneA:a5lronly:+"]
    assert _read(str(out / "LR.A3SS.lst")) == ["chr2:geneB:a3one:-"]
    assert _read(str(out / "LR.IR.lst")) == ["chr3:geneC:ir1:+"]

    # confident.* are relative symlinks to the LR.* lists (no NGS -> same content)
    assert os.path.islink(str(out / "confident.SE.lst"))
    assert os.readlink(str(out / "confident.SE.lst")) == "LR.SE.lst"
    assert _read(str(out / "confident.SE.lst")) == ["chr1:geneA:se1:+"]
    assert _read(str(out / "confident.A5SS.lst")) == _read(str(out / "LR.A5SS.lst"))
    assert _read(str(out / "confident.A3SS.lst")) == _read(str(out / "LR.A3SS.lst"))

    # no NGS / APA artifacts produced
    assert not (out / "NGS.SE.lst").exists()
    assert not (out / "Common.SE.lst").exists()
    assert not (out / "LR.APA.lst").exists()


def test_build_anno_as_lists_with_pa(tmp_path):
    ordinary = _write_ordinary_as_lr(tmp_path)
    pa = tmp_path / "reads.paGrouped.tsv"
    pa.write_text("chr3\t+\tgeneC\t100,200\nchr3\t+\tgeneD\t500\n")  # geneD single PA -> skip
    out = tmp_path / "char_pa"
    build_anno_as_lists(str(ordinary), str(out), pa_grouped=str(pa), ngs=False)
    assert _read(str(out / "LR.APA.lst")) == ["chr3:+:geneC:100,200"]


def test_build_anno_as_lists_ngs_branch_and_venn_contract(tmp_path):
    ordinary = _write_ordinary_as_lr(tmp_path)
    ngsDir = tmp_path / "ordinary_as" / "NGS"
    ngsDir.mkdir(parents=True)
    (ngsDir / "SE.bed12+").write_text(
        "chr1\t1000\t2000\tgeneA:se1\t.\t+\tx\n"        # shared with LR
        "chr1\t1000\t2000\tgeneA:seNGSonly\t.\t+\tx\n"  # NGS-only
    )
    (ngsDir / "A5SS.known.bed6+").write_text("chr1\t100\t200\tgeneA:a5known\t.\t+\tg\n")
    (ngsDir / "A5SS.LR.bed6+").write_text("chr1\t100\t200\tgeneA:a5shared\t.\t+\tg\n")
    (ngsDir / "A3SS.known.bed6+").write_text("chr2\t500\t600\tgeneB:a3known\t.\t-\tg\n")
    (ngsDir / "A3SS.LR.bed6+").write_text("chr2\t500\t600\tgeneB:a3one\t.\t-\tg\n")

    runner = _RecordingRunner()
    out = tmp_path / "char_ngs"
    build_anno_as_lists(str(ordinary), str(out), ngs=True, runner=runner)

    # BUG-117 fixed: known events remain when LR events are appended.
    assert _read(str(out / "NGS.A5SS.lst")) == [
        "chr1:geneA:a5known:+",
        "chr1:geneA:a5shared:+",
    ]
    assert _read(str(out / "NGS.A3SS.lst")) == [
        "chr2:geneB:a3known:-",
        "chr2:geneB:a3one:-",
    ]

    # Common.* = LR rows whose event is also called by NGS
    assert _read(str(out / "Common.SE.lst")) == ["chr1:geneA:se1:+"]
    assert _read(str(out / "Common.A5SS.lst")) == ["chr1:geneA:a5shared:+"]

    # confident.* now point at the Common.* intersections
    assert os.path.islink(str(out / "confident.SE.lst"))
    assert os.readlink(str(out / "confident.SE.lst")) == "Common.SE.lst"
    assert _read(str(out / "confident.SE.lst")) == ["chr1:geneA:se1:+"]

    # LR.specific.* = LR events NGS did not call
    assert _read(str(out / "LR.specific.SE.lst")) == []
    assert _read(str(out / "LR.specific.A5SS.lst")) == ["chr1:geneA:a5lronly:+"]

    # venn.R command contract: one NGS-vs-LR overlap plot per AS type, run by Rscript in the R env, non-fatal
    assert len(runner.calls) == 3
    for (cmd, kwargs), asType in zip(runner.calls, ("SE", "A5SS", "A3SS")):
        assert cmd.startswith("Rscript ") and cmd.split()[1].endswith("/venn.R")
        assert cmd.split()[2:4] == ["NGS.{}.lst".format(asType), "LR.{}.lst".format(asType)]
        assert "-p={}.pdf".format(asType) in cmd
        assert kwargs["conda_env"] == "iflas-diffsplice" and kwargs["check"] is False and kwargs["cwd"] == str(out)


# ---------------------------------------------------------------------------
# characterize_sample  (file-level charaterize_as driver)
# ---------------------------------------------------------------------------

def _write_base_ordinary_as(base):
    """Write a minimal <base>/as_events/ordinary_as/LR tree (find_as contract)."""
    lrDir = base / "as_events" / "ordinary_as" / "LR"
    lrDir.mkdir(parents=True)
    (lrDir / "SE.bed12+").write_text("chr1\t1000\t2000\tgeneA:se1\t.\t+\textra\n")
    (lrDir / "A5SS.bed6+").write_text("chr1\t100\t200\tgeneA:a5one\t.\t+\tgeneA\n")
    (lrDir / "A3SS.bed6+").write_text("chr2\t500\t600\tgeneB:a3one\t.\t-\tgeneB\n")
    # IR event id is Gene:RetainedIntronStart-RetainedIntronEnd
    (lrDir / "IR.bed6+").write_text("chr3\t100\t200\tgeneC:100-200\t.\t+\tgeneC\n")
    return base


def test_characterize_sample_lists_and_splicesite(tmp_path):
    base = _write_base_ordinary_as(tmp_path)
    runner = _RecordingRunner()
    prev = os.getcwd()
    out = characterize_sample(str(base), runner=runner)

    # working directory restored after the driver finishes
    assert os.getcwd() == prev

    charDir = tmp_path / "as_events" / "characterization"
    assert out["out_dir"] == str(charDir)

    # front-half lists + confident links produced via build_anno_as_lists
    assert _read(str(charDir / "LR.SE.lst")) == ["chr1:geneA:se1:+"]
    assert os.path.islink(str(charDir / "confident.SE.lst"))
    assert os.readlink(str(charDir / "confident.SE.lst")) == "LR.SE.lst"

    # IR splice-site is the pure branch -> real deterministic output
    assert _read(str(charDir / "IR.splicesite")) == ["chr3\t100\t200\t+"]

    # SE/A3SS/A5SS splice-site decompose routed through the runner (Perl helpers)
    perl_calls = [c for (c, _) in runner.calls if c and c[0] == "perl"]
    assert len(perl_calls) == 3
    assert os.path.basename(perl_calls[0][1]) == "seDecompose.pl"
    assert perl_calls[0][-1] == "confident.SE.lst"
    assert os.path.basename(perl_calls[1][1]) == "anssDecompose.pl"
    assert perl_calls[1][2:4] == ["-n", "3"]
    assert perl_calls[1][-1] == "confident.A3SS.lst"
    assert perl_calls[2][2:4] == ["-n", "5"]
    assert perl_calls[2][-1] == "confident.A5SS.lst"


def test_characterize_sample_missing_ordinary_as(tmp_path):
    with pytest.raises(FileNotFoundError):
        characterize_sample(str(tmp_path))


# ---------------------------------------------------------------------------
# back half: splice-site dinucleotide profiles
# ---------------------------------------------------------------------------

@pytest.mark.skipif(shutil.which("bedtools") is None, reason="needs bedtools (pybedtools getfasta)")
def test_splicesite2seq_counts_dinucleotide_pairs(tmp_path):
    # intron 1: GT...AG on +; intron 2: GC...AG on +; intron 3 on -: genomic CT...AC reads as GT-AG
    genome = tmp_path / "genome.fa"
    seq = list("A" * 200)
    seq[50:52], seq[78:80] = "GT", "AG"
    seq[100:102], seq[128:130] = "GC", "AG"
    seq[150:152], seq[178:180] = "CT", "AC"
    genome.write_text(">chr1\n" + "".join(seq) + "\n")
    sites = tmp_path / "x.splicesite"
    sites.write_text("chr1\t50\t80\t+\nchr1\t100\t130\t+\nchr1\t150\t180\t-\n")
    splicesite2seq(str(genome), str(sites), outFile=str(tmp_path / "x.tmp"))
    assert _read(str(tmp_path / "x.tmp")) == ["GT-AG\t2", "GC-AG\t1"]
    splicesite2seq(str(genome), str(sites), noThreeEnd=True, outFile=str(tmp_path / "donor.tmp"))
    assert _read(str(tmp_path / "donor.tmp")) == ["GT\t2", "GC\t1"]


def test_draw_ss_motif_normalises_top3_plus_other(tmp_path, monkeypatch):
    motif = tmp_path / "SE.inc.tmp"
    motif.write_text("GT-AG\t6\nGC-AG\t2\nAT-AC\t1\nGA-AG\t1\n")
    runner = _RecordingRunner()
    monkeypatch.chdir(tmp_path)
    drawSSmotif(asMotif=str(motif), outPrefix="SE.inc", asType="A5SS", runner=runner)
    assert _read(str(tmp_path / "SE.inc.tmp.txt")) == ["GT-AG-\t0.6", "GC-AG-\t0.2", "AT-AC-\t0.1", "Other\t0.1"]
    cmd = runner.calls[-1][0]
    assert "bar.R" in str(cmd) and "-p=SE.inc.ssMotif.pdf" in str(cmd)
