import subprocess
from pathlib import Path

from iflas.cli import main
from pypdf import PdfReader, PdfWriter


def _write_report_config(tmp_path: Path, out_dir: Path) -> tuple[Path, Path]:
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        f"""\
reference:
  - ref_strain: B73
    ref_genome: {tmp_path / "genome.fa"}
    ref_gtf: {tmp_path / "genes.gtf"}
dirs:
  out_dir: {out_dir}
"""
    )
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample_name\tproject_name\tref_strain\n"
        "s1\tdemo\tB73\n"
    )
    return cfg, samples


def test_iflas_report_all_expands_to_migrated_report_helpers(tmp_path, monkeypatch, capsys):
    import iflas.report.html as report_html
    import iflas.report.plots as report_plots

    out_dir = tmp_path / "out"
    cfg, samples = _write_report_config(tmp_path, out_dir)
    calls = []

    def fake_sample_report(name):
        def _fake(config, *, sample_names=None, runner=None):
            calls.append((name, list(sample_names or []), runner))
            return {"plots": [f"{name}.pdf"], "by_sample": {}}

        return _fake

    monkeypatch.setattr(report_plots, "plot_basic_report", fake_sample_report("basic"))
    monkeypatch.setattr(report_plots, "plot_as_pattern_report", fake_sample_report("asp"))
    monkeypatch.setattr(report_plots, "plot_target_gene_structure_report", fake_sample_report("geneStruc"))
    monkeypatch.setattr(report_plots, "plot_novel_hq_as_report", fake_sample_report("novelHqAs"))
    monkeypatch.setattr(report_plots, "plot_asas_report", fake_sample_report("asas"))
    monkeypatch.setattr(report_plots, "plot_palen_report", fake_sample_report("palen"))

    def fake_diff_report(config, *, runner=None):
        calls.append(("diff", [], runner))
        return {"plots": ["diff.pdf"], "by_comparison": {}}

    monkeypatch.setattr(report_plots, "plot_diff_as_report", fake_diff_report)

    def fake_html_report(config, *, sample_names=None):
        calls.append(("html", list(sample_names or []), None))
        return {"main_page": str(out_dir / "reports" / "html" / "iflas_report.html"), "json": str(out_dir / "reports" / "json" / "search.json")}

    monkeypatch.setattr(report_html, "generate_html_report", fake_html_report)

    rc = main(["report", "-all", "-c", str(cfg), "-s", str(samples), "--sample", "s1"])

    assert rc == 0
    assert calls == [
        ("basic", ["s1"], None),
        ("asp", ["s1"], None),
        ("geneStruc", ["s1"], None),
        ("novelHqAs", ["s1"], None),
        ("diff", [], None),
        ("asas", ["s1"], None),
        ("palen", ["s1"], None),
        ("html", ["s1"], None),
    ]
    out = capsys.readouterr().out
    assert "Generated basic report plots" in out
    assert "Generated AS pattern report plots" in out
    assert "Generated target gene structure report plots" in out
    assert "Generated differential AS report plots" in out
    assert "Generated allele-specific AS report plots" in out
    assert "Generated poly(A) tail report plots" in out
    assert "Wrote HTML report" in out


def _bed12plus(name, exons, group, chrom="chr1", strand="+", gene="geneX"):
    start = exons[0][0]
    end = exons[-1][1]
    sizes = ",".join(str(e - s) for s, e in exons)
    rel_starts = ",".join(str(s - start) for s, e in exons)
    fields = [
        chrom,
        start,
        end,
        name,
        0,
        strand,
        start,
        end,
        "0,0,0",
        len(exons),
        sizes,
        rel_starts,
        gene,
        gene,
        group,
    ]
    return "\t".join(map(str, fields))


def test_iflas_report_novel_hq_as_routes_plot_through_rscript(tmp_path, monkeypatch, capsys):
    import iflas.report.plots as report_plots

    out_dir = tmp_path / "out"
    cfg, samples = _write_report_config(tmp_path, out_dir)
    score = out_dir / "demo" / "s1" / "isoIdentity" / "isoformScore.txt"
    score.parent.mkdir(parents=True)
    score.write_text("geneA\tPB.1.1\t2\t4\t0.5\tNovel\n")

    calls = []

    def fake_run_r_function(source, function_name, args, cwd=None, runner=None):
        calls.append((source, function_name, args, cwd, runner))
        Path(cwd, args[1]).write_text("fake pdf")
        return subprocess.CompletedProcess(["Rscript"], 0)

    monkeypatch.setattr(report_plots, "run_r_function", fake_run_r_function)

    rc = main(["report", "--novel-hq-as", "-c", str(cfg), "-s", str(samples)])

    assert rc == 0
    assert len(calls) == 1
    source, function_name, r_args, cwd, runner = calls[0]
    assert "plotNovelHqAS <- function" in source
    assert function_name == "plotNovelHqAS"
    assert r_args == [str(score.resolve()), "isoformScore.pdf"]
    assert cwd == str((out_dir / "reports" / "demo_s1").resolve())
    assert runner is None
    assert (out_dir / "reports" / "demo_s1" / "isoformScore.pdf").is_file()
    assert "Generated report plots" in capsys.readouterr().out


def test_iflas_report_diff_routes_summary_plot_through_rscript(tmp_path, monkeypatch, capsys):
    import iflas.report.plots as report_plots

    out_dir = tmp_path / "out"
    cfg, samples = _write_report_config(tmp_path, out_dir)
    sig_dir = out_dir / "das" / "condA_vs_condB.sigDiffAS"
    sig_dir.mkdir(parents=True)
    (sig_dir / "IR.sig.txt").write_text("event\tpvalue\nir1\t0.01\n")
    (sig_dir / "SE.sig.txt").write_text("event\tpvalue\nse1\t0.01\nse2\t0.02\n")
    (sig_dir / "A3SS.sig.txt").write_text("event\tpvalue\n")
    (sig_dir / "A5SS.sig.txt").write_text("event\tpvalue\na5ss1\t0.03\n")

    calls = []

    def fake_run_r_function(source, function_name, args, cwd=None, runner=None):
        calls.append((source, function_name, args, cwd, runner))
        Path(cwd, args[1]).write_text("fake pdf")
        return subprocess.CompletedProcess(["Rscript"], 0)

    monkeypatch.setattr(report_plots, "run_r_function", fake_run_r_function)

    rc = main(["report", "-diff", "-c", str(cfg), "-s", str(samples)])

    assert rc == 0
    summary = sig_dir / "sigDiff.AS.txt"
    assert summary.read_text() == "AS_type\tCount\nIR\t1\nSE\t2\nA3SS\t0\nA5SS\t1\n"
    assert len(calls) == 1
    source, function_name, r_args, cwd, runner = calls[0]
    assert "plotDiffAS <- function" in source
    assert function_name == "plotDiffAS"
    assert r_args == ["sigDiff.AS.txt", "sigDiff.AS_distribution.pdf"]
    assert cwd == str(sig_dir.resolve())
    assert runner is None
    assert (sig_dir / "sigDiff.AS_distribution.pdf").is_file()
    assert "Generated differential AS report plots" in capsys.readouterr().out


def test_iflas_report_basic_routes_legacy_plots(tmp_path, monkeypatch, capsys):
    import iflas.report.plots as report_plots

    out_dir = tmp_path / "out"
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\n" + "A" * 1000 + "\n")
    gtf = tmp_path / "genes.gtf"
    gtf.write_text("")
    ref_gpe = tmp_path / "genes.gpe"
    ref_gpe.write_text("tx1\tchr1\t+\t0\t100\t0\t100\t1\t0,\t100,\t0\tgeneA\tcmpl\tcmpl\t0,\n")
    flnc = tmp_path / "flnc.fa"
    flnc.write_text(">read1\nACGTACGT\n>read2\nGGGGAAAA\n")

    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        f"""\
reference:
  - ref_strain: B73
    ref_genome: {genome}
    ref_gtf: {gtf}
    ref_gpe: {ref_gpe}
dirs:
  out_dir: {out_dir}
"""
    )
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample_name\tproject_name\tref_strain\ttgs_plat\tdata_processed_location\tuse_fmlrc2\n"
        f"s1\tdemo\tB73\tpacbio\t{flnc}\ttrue\n"
    )

    mapping_dir = out_dir / "demo" / "s1" / "mapping"
    mapping_dir.mkdir(parents=True)
    (mapping_dir / "rawFlnc.addCVandID.bed12+").write_text("chr1\t0\t10\tread1\t0\t+\t0\t10\t0\t1\t10,\t0,\t90\t1\n")
    (mapping_dir / "flnc.addCVandID.bed12+").write_text("chr1\t0\t10\tread1\t0\t+\t0\t10\t0\t1\t10,\t0,\t95\t1\n")

    calls = []

    def fake_run_r_function(source, function_name, args, cwd=None, runner=None):
        calls.append(("r", source, function_name, args, cwd, runner))
        Path(cwd, args[-1]).write_text("fake pdf")
        return subprocess.CompletedProcess(["Rscript"], 0)

    class DummyRunner:
        def __init__(self):
            self.commands = []

        def run(self, cmd, **kwargs):
            self.commands.append((cmd, kwargs))
            if cmd[:2] == ["seqkit", "head"]:
                return subprocess.CompletedProcess(cmd, 0, stdout=flnc.read_text(), stderr="")
            if cmd[:2] == ["seqkit", "stat"]:
                return subprocess.CompletedProcess(cmd, 0, stdout="file\tformat\nflnc.fa\tFASTA\n", stderr="")

            script = cmd[-1] if isinstance(cmd, list) and cmd[:2] == ["bash", "-c"] else str(cmd)
            cwd = Path(kwargs["cwd"])
            if "GC_of_raw_flnc_reads.log" in script and "seqkit fx2tab -n -g" in script:
                (cwd / "GC_of_raw_flnc_reads.log").write_text("read1\t50\nread2\t50\n")
            if "-p=GC_of_raw_flnc_reads.pdf" in script:
                (cwd / "GC_of_raw_flnc_reads.pdf").write_text("fake pdf")
            if "-p=GC_across_raw_flnc_read.pdf" in script:
                (cwd / "GC_across_raw_flnc_read.pdf").write_text("fake pdf")
            if "readsLength.lst" in script and "seqkit fx2tab -n -l" in script:
                (cwd / "readsLength.lst").write_text("8\n8\n")
            if "refGeneLength.lst" in script:
                (cwd / "refGeneLength.lst").write_text("100\n")
            if "-p=LengthDistribution.curve.pdf" in script:
                (cwd / "LengthDistribution.curve.pdf").write_text("fake pdf")
            if "-p=LengthDistribution.box.pdf" in script:
                (cwd / "LengthDistribution.box.pdf").write_text("fake pdf")
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    dummy_runner = DummyRunner()
    monkeypatch.setattr(report_plots, "run_r_function", fake_run_r_function)
    monkeypatch.setattr(report_plots, "ToolRunner", lambda: dummy_runner)

    rc = main(["report", "-basic", "-c", str(cfg), "-s", str(samples)])

    assert rc == 0
    basic_dir = out_dir / "reports" / "demo_s1" / "basicStatistics"
    for name in (
        "readsCorrectResult.pdf",
        "GC_of_raw_flnc_reads.pdf",
        "GC_across_raw_flnc_read.pdf",
        "LengthDistribution.curve.pdf",
        "LengthDistribution.box.pdf",
    ):
        assert (basic_dir / name).is_file()
    assert calls[0][2] == "plotReadsCorrectedEval"
    assert calls[0][3] == [
        str((mapping_dir / "rawFlnc.addCVandID.bed12+").resolve()),
        str((mapping_dir / "flnc.addCVandID.bed12+").resolve()),
        "readsCorrectResult.pdf",
    ]
    scripts = "\n".join(cmd[-1] for cmd, _ in dummy_runner.commands if isinstance(cmd, list) and cmd[:2] == ["bash", "-c"])
    assert "seqkit fx2tab -n -g" in scripts
    assert "distrCurve.R" in scripts
    assert "gpe2bed.pl" in scripts
    assert "LengthDistribution.curve.pdf" in scripts
    assert "Generated basic report plots" in capsys.readouterr().out


def test_gc_across_read_keeps_legacy_integer_division(tmp_path, monkeypatch):
    import iflas.report.plots as report_plots

    flnc = tmp_path / "reads.fa"
    out = tmp_path / "gc.log"
    flnc.write_text(">read1\n" + "G" * 41 + "\n")
    monkeypatch.setattr(report_plots, "validateFaAndFqFile", lambda *args, **kwargs: "fasta")

    report_plots._write_gc_across_read(flnc, out, interval=20, runner=object())

    fields = out.read_text().strip().split("\t")
    assert fields[0] == "read1"
    assert len(fields[1:]) == 21


def test_iflas_report_asp_routes_as_pattern_plots(tmp_path, monkeypatch, capsys):
    import iflas.report.plots as report_plots

    out_dir = tmp_path / "out"
    cfg, samples = _write_report_config(tmp_path, out_dir)
    charac_dir = out_dir / "demo" / "s1" / "as_events" / "characterization"
    charac_dir.mkdir(parents=True)
    for name, rows in {
        "IR.anno.lst": ["ir_anno_1"],
        "IR.novel.lst": ["ir_novel_1", "ir_novel_2"],
        "SE.anno.lst": ["se_anno_1"],
        "SE.novel.lst": ["se_novel_1", "se_novel_2", "se_novel_3"],
        "A3SS.anno.lst": ["a3_anno_1"],
        "A3SS.novel.lst": ["a3_novel_1"],
        "A5SS.anno.lst": ["a5_anno_1", "a5_anno_2"],
        "A5SS.novel.lst": ["a5_novel_1"],
        "annoPA.tsv": ["pa_anno_1"],
        "novelPA.tsv": ["pa_novel_1", "pa_novel_2"],
    }.items():
        (charac_dir / name).write_text("\n".join(rows) + "\n")
    (charac_dir / "statistic.APA.tsv").write_text(
        "gene1\tchr1_100\tKnown\n"
        "gene2\tchr1_200\tNovel\n"
        "gene3\tchr1_300\tNovel\n"
    )
    for name, rows in {
        "IR.tmp.txt": ["GT-AG\t2"],
        "SE.inc.tmp.txt": ["GT-AG\t3"],
        "SE.exc.tmp.txt": ["GC-AG\t1"],
        "A3SS.inc.tmp.txt": ["GT-AG\t4"],
        "A3SS.exc.tmp.txt": ["AT-AC\t1"],
        "A5SS.inc.tmp.txt": ["GT-AG\t5"],
        "A5SS.exc.tmp.txt": ["GC-AG\t2"],
    }.items():
        (charac_dir / name).write_text("\n".join(rows) + "\n")

    calls = []

    def fake_run_r_function(source, function_name, args, cwd=None, runner=None):
        calls.append((source, function_name, args, cwd, runner))
        Path(cwd, args[1]).write_text("fake pdf")
        return subprocess.CompletedProcess(["Rscript"], 0)

    monkeypatch.setattr(report_plots, "run_r_function", fake_run_r_function)

    rc = main(["report", "-asp", "-c", str(cfg), "-s", str(samples)])

    assert rc == 0
    as_pattern_dir = out_dir / "reports" / "demo_s1" / "asPattern"
    anno_file = as_pattern_dir / "demo_s1.AS_annotation.txt"
    splice_file = as_pattern_dir / "demo_s1.AS_spliceSite.txt"
    assert anno_file.read_text() == (
        "AS_type\tAnnotation\tCount\n"
        "IR\tAnno\t1\n"
        "IR\tNovel\t2\n"
        "SE\tAnno\t1\n"
        "SE\tNovel\t3\n"
        "A3SS\tAnno\t1\n"
        "A3SS\tNovel\t1\n"
        "A5SS\tAnno\t2\n"
        "A5SS\tNovel\t1\n"
        "APA\tAnno\t1\n"
        "APA\tNovel\t2\n"
        "PA\tAnno\t1\n"
        "PA\tNovel\t2\n"
    )
    assert splice_file.read_text() == (
        "Dinucleotide\tFrequency\tAS_type\tCategory\n"
        "GT-AG\t2\tIR\tInc\n"
        "GT-AG\t3\tSE\tInc\n"
        "GC-AG\t1\tSE\tExc\n"
        "GT-AG\t4\tA3SS\tInc\n"
        "AT-AC\t1\tA3SS\tExc\n"
        "GT-AG\t5\tA5SS\tInc\n"
        "GC-AG\t2\tA5SS\tExc\n"
    )
    assert len(calls) == 2
    source, function_name, r_args, cwd, runner = calls[0]
    assert "plotAsCountStatistics <- function" in source
    assert function_name == "plotAsCountStatistics"
    assert r_args == ["demo_s1.AS_annotation.txt", "demo_s1.AS_annotation.pdf"]
    assert cwd == str(as_pattern_dir.resolve())
    assert runner is None
    source, function_name, r_args, cwd, runner = calls[1]
    assert "plotAsDinucleotideStatistics <- function" in source
    assert function_name == "plotAsDinucleotideStatistics"
    assert r_args == ["demo_s1.AS_spliceSite.txt", "demo_s1.AS_spliceSite.pdf"]
    assert cwd == str(as_pattern_dir.resolve())
    assert runner is None
    assert (as_pattern_dir / "demo_s1.AS_annotation.pdf").is_file()
    assert (as_pattern_dir / "demo_s1.AS_spliceSite.pdf").is_file()
    assert "Generated AS pattern report plots" in capsys.readouterr().out


def test_iflas_report_asp_tolerates_missing_pa_and_splice_site_back_half(tmp_path, monkeypatch, capsys):
    import iflas.report.plots as report_plots

    out_dir = tmp_path / "out"
    cfg, samples = _write_report_config(tmp_path, out_dir)
    charac_dir = out_dir / "demo" / "s1" / "as_events" / "characterization"
    charac_dir.mkdir(parents=True)
    for name, rows in {
        "IR.anno.lst": ["ir_anno_1"],
        "IR.novel.lst": [],
        "SE.anno.lst": ["se_anno_1"],
        "SE.novel.lst": ["se_novel_1"],
        "A3SS.anno.lst": [],
        "A3SS.novel.lst": [],
        "A5SS.anno.lst": [],
        "A5SS.novel.lst": ["a5_novel_1"],
    }.items():
        (charac_dir / name).write_text("\n".join(rows) + ("\n" if rows else ""))

    calls = []

    def fake_run_r_function(source, function_name, args, cwd=None, runner=None):
        calls.append((source, function_name, args, cwd, runner))
        Path(cwd, args[1]).write_text("fake pdf")
        return subprocess.CompletedProcess(["Rscript"], 0)

    monkeypatch.setattr(report_plots, "run_r_function", fake_run_r_function)

    rc = main(["report", "-asp", "-c", str(cfg), "-s", str(samples)])

    assert rc == 0
    as_pattern_dir = out_dir / "reports" / "demo_s1" / "asPattern"
    anno_file = as_pattern_dir / "demo_s1.AS_annotation.txt"
    assert anno_file.read_text() == (
        "AS_type\tAnnotation\tCount\n"
        "IR\tAnno\t1\n"
        "IR\tNovel\t0\n"
        "SE\tAnno\t1\n"
        "SE\tNovel\t1\n"
        "A3SS\tAnno\t0\n"
        "A3SS\tNovel\t0\n"
        "A5SS\tAnno\t0\n"
        "A5SS\tNovel\t1\n"
        "APA\tAnno\t0\n"
        "APA\tNovel\t0\n"
        "PA\tAnno\t0\n"
        "PA\tNovel\t0\n"
    )
    assert len(calls) == 1
    assert calls[0][1] == "plotAsCountStatistics"
    assert (as_pattern_dir / "demo_s1.AS_annotation.pdf").is_file()
    assert not (as_pattern_dir / "demo_s1.AS_spliceSite.txt").exists()
    assert "Generated AS pattern report plots" in capsys.readouterr().out


def test_iflas_report_gene_struc_merges_existing_visual_as_pdfs(tmp_path, capsys):
    out_dir = tmp_path / "out"
    cfg, samples = _write_report_config(tmp_path, out_dir)

    def write_pdf(path: Path, pages: int) -> None:
        path.parent.mkdir(parents=True)
        writer = PdfWriter()
        for _ in range(pages):
            writer.add_blank_page(width=72, height=72)
        with path.open("wb") as handle:
            writer.write(handle)

    write_pdf(out_dir / "demo" / "s1" / "visual_as" / "geneA" / "geneA.pdf", 1)
    write_pdf(out_dir / "demo" / "s1" / "isoViewer" / "geneB" / "geneB.pdf", 2)

    rc = main(["report", "-geneStruc", "-c", str(cfg), "-s", str(samples)])

    assert rc == 0
    merged = out_dir / "reports" / "demo_s1" / "geneStrucPlots" / "allTargetGeneStructure.pdf"
    assert merged.is_file()
    with merged.open("rb") as handle:
        assert len(PdfReader(handle).pages) == 3
    assert "Generated target gene structure report plots" in capsys.readouterr().out


def test_iflas_report_palen_routes_as_and_apa_plots(tmp_path, monkeypatch, capsys):
    import iflas.report.plots as report_plots

    out_dir = tmp_path / "out"
    gtf = tmp_path / "genes.gtf"
    gtf.write_text("")
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        f"""\
reference:
  - ref_strain: B73
    ref_genome: {tmp_path / "genome.fa"}
    ref_gtf: {gtf}
dirs:
  out_dir: {out_dir}
"""
    )
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample_name\tproject_name\tref_strain\n"
        "s1\tdemo\tB73\n"
    )

    base = out_dir / "demo" / "s1"
    merge_dir = base / "palenAS" / "mergeByJunc"
    apa_dir = base / "palenAS" / "palenAPA"
    refine_dir = base / "refine"
    mapping_dir = base / "mapping"
    for path in (merge_dir, apa_dir, refine_dir, mapping_dir):
        path.mkdir(parents=True)

    sig_line = "\t".join(
        [
            "chr1",
            "100",
            "200",
            "PB.1.1",
            "0",
            "+",
            "100",
            "200",
            "0,0,0",
            "1",
            "100",
            "0",
            "geneX",
            "event1",
            "inclusion",
            "10,11",
            "2",
        ]
    )
    for as_type in ("IR", "SE", "A3SS", "A5SS"):
        (merge_dir / f"{as_type}.palenAndAS.sig.bed12+").write_text("")
    (merge_dir / "SE.palenAndAS.sig.bed12+").write_text(sig_line + "\n" + sig_line + "\n")

    apa_line = "\t".join(
        [
            "geneX",
            "+",
            "chr1_110;chr1_210",
            "6;6",
            "10,11,12,13,14,15;50,51,52,53,54,55",
        ]
    )
    (apa_dir / "apaRelatedPalen.txt").write_text(apa_line + "\n")
    (refine_dir / "reads.assigned.unambi.bed12+").write_text(
        "chr1\t100\t400\tPB.1.1\t0\t+\t100\t400\t0,0,0\t2\t100,100\t0,200\tPB.1.1\tNA\tgeneX\n"
    )
    align_bam = mapping_dir / "flnc.mm2.sorted.bam"
    align_bam.write_text("")

    calls = []

    def fake_run_r_function(source, function_name, args, cwd=None, runner=None):
        calls.append((source, function_name, args, cwd, runner))
        if function_name == "plotPaTailAsStructure":
            Path(cwd, "geneX.palenAndAS.pdf").write_text("fake as pdf")
        else:
            Path(cwd, f"{args[0]}.palenAndAPA.pdf").write_text("fake apa pdf")
        return subprocess.CompletedProcess(["Rscript"], 0)

    class DummyRunner:
        def __init__(self):
            self.commands = []

        def run(self, cmd, **kwargs):
            self.commands.append((cmd, kwargs))
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    dummy_runner = DummyRunner()
    monkeypatch.setattr(report_plots, "run_r_function", fake_run_r_function)
    monkeypatch.setattr(report_plots, "run_r_calls", lambda source, name, calls_, cwd=None, runner=None:
                        [fake_run_r_function(source, name, args, cwd=cwd, runner=runner) for args in calls_])
    monkeypatch.setattr(report_plots, "ToolRunner", lambda: dummy_runner)

    rc = main(["report", "-palen", "-c", str(cfg), "-s", str(samples)])

    assert rc == 0
    report_dir = out_dir / "reports" / "demo_s1"
    merged = report_dir / "palenAsPlots" / "palenAndAS.sig.bed12+"
    assert merged.read_text() == sig_line + "\n"
    assert len(calls) == 2
    source, function_name, r_args, cwd, runner = calls[0]
    assert "plotPaTailAsStructure <- function" in source
    assert function_name == "plotPaTailAsStructure"
    assert r_args == ["palenAndAS.sig.bed12+"]
    assert cwd == str((report_dir / "palenAsPlots").resolve())
    assert runner is dummy_runner
    source, function_name, r_args, cwd, runner = calls[1]
    assert "plotPaTailApaStructure <- function" in source
    assert "makeTxDbFromGFF" in source
    assert str(gtf) in source
    assert function_name == "plotPaTailApaStructure"
    assert r_args == [
        "geneX",
        "chr1_110;chr1_210",
        "10,11,12,13,14,15;50,51,52,53,54,55",
        "6;6",
        str(align_bam.resolve()),
        "chr1",
        100,
        400,
    ]
    assert cwd == str((report_dir / "palenApaPlots").resolve())
    assert runner is dummy_runner
    assert dummy_runner.commands == [(["samtools", "index", str(align_bam.resolve())], {})]
    assert (report_dir / "palenAsPlots" / "geneX.palenAndAS.pdf").is_file()
    assert (report_dir / "palenApaPlots" / "geneX.palenAndAPA.pdf").is_file()
    assert "Generated poly(A) tail report plots" in capsys.readouterr().out


def test_iflas_report_asas_routes_allele_as_plots(tmp_path, monkeypatch, capsys):
    import iflas.report.plots as report_plots

    out_dir = tmp_path / "out"
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\n" + "A" * 1000 + "\n")
    cfg, samples = _write_report_config(tmp_path, out_dir)

    base = out_dir / "demo" / "s1"
    for path in (
        base / "asas",
        base / "hqIsoforms",
        base / "collapse",
        base / "mapping",
    ):
        path.mkdir(parents=True)
    (base / "asas" / "partialAsRelatedHaplotype.txt").write_text(
        "\t".join(["PB.1", "SE", "PB.1", "PB.1:150@200-250@300", "hapREF", "PB.1.1", "hapALT", "PB.1.2"]) + "\n"
    )
    (base / "hqIsoforms" / "hq.collapsed.bed12+").write_text(
        "\n".join(
            [
                _bed12plus("PB.1.1", [(100, 150), (200, 250), (300, 350)], "PB.1.1"),
                _bed12plus("PB.1.2", [(100, 150), (300, 350)], "PB.1.2"),
            ]
        )
        + "\n"
    )
    (base / "collapse" / "tofu.collapsed.group.txt").write_text(
        "PB.1.1\tread_inc1,read_inc2\n"
        "PB.1.2\tread_exc1,read_exc2\n"
    )
    flnc_bam = base / "mapping" / "flnc.mm2.sorted.bam"
    flnc_bam.write_text("")

    calls = []

    def fake_run_r_function(source, function_name, args, cwd=None, runner=None):
        calls.append((source, function_name, args, cwd, runner))
        writer = PdfWriter()
        writer.add_blank_page(width=72, height=72)
        with open(f"{args[-1]}.pdf", "wb") as handle:
            writer.write(handle)
        return subprocess.CompletedProcess(["Rscript"], 0)

    class DummyRunner:
        def __init__(self):
            self.commands = []

        def run(self, cmd, **kwargs):
            self.commands.append((cmd, kwargs))
            script = cmd[-1] if isinstance(cmd, list) and cmd[:2] == ["bash", "-c"] else ""
            cwd = Path(kwargs.get("cwd", "."))
            if "haplo1isosOut.gpe" in script:
                (cwd / "haplo1isosOut.gpe").write_text("gpe1\n")
            if "haplo2isosOut.gpe" in script:
                (cwd / "haplo2isosOut.gpe").write_text("gpe2\n")
            if "haplo1isosOut.gtf" in script:
                (cwd / "haplo1isosOut.gtf").write_text("gtf1\n")
            if "haplo2isosOut.gtf" in script:
                (cwd / "haplo2isosOut.gtf").write_text("gtf2\n")
            if "haplo1.flnc.sorted.bam" in script:
                (cwd / "haplo1.flnc.sorted.bam").write_text("")
            if "haplo2.flnc.sorted.bam" in script:
                (cwd / "haplo2.flnc.sorted.bam").write_text("")
            if "PB.1.PB.1.SE.allele_as.flnc.sorted.bam" in script:
                (cwd / "PB.1.PB.1.SE.allele_as.flnc.sorted.bam").write_text("")
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    dummy_runner = DummyRunner()
    monkeypatch.setattr(report_plots, "run_r_function", fake_run_r_function)
    monkeypatch.setattr(report_plots, "run_r_calls", lambda source, name, calls_, cwd=None, runner=None:
                        [fake_run_r_function(source, name, args, cwd=cwd, runner=runner) for args in calls_])
    monkeypatch.setattr(report_plots, "ToolRunner", lambda: dummy_runner)

    rc = main(["report", "-asas", "-c", str(cfg), "-s", str(samples)])

    assert rc == 0
    report_dir = out_dir / "reports" / "demo_s1" / "alleleAsPlots"
    event_pdf = report_dir / "PB.1.PB.1.SE.allele_as.pdf"
    assert event_pdf.is_file()
    assert (out_dir / "reports" / "demo_s1" / "alleleAS.pdf").is_file()
    assert len(calls) == 1
    source, function_name, r_args, cwd, runner = calls[0]
    assert "plotAlleleAsStructure <- function" in source
    assert function_name == "plotAlleleAsStructure"
    assert r_args == [
        str(genome.resolve()),
        f"{report_dir / 'PB.1.PB.1.SE.allele_as' / 'haplo1isosOut.gtf'},{report_dir / 'PB.1.PB.1.SE.allele_as' / 'haplo2isosOut.gtf'}",
        str(report_dir / "PB.1.PB.1.SE.allele_as" / "PB.1.PB.1.SE.allele_as.flnc.sorted.bam"),
        f"{report_dir / 'PB.1.PB.1.SE.allele_as' / 'haplo1.flnc.sorted.bam'},{report_dir / 'PB.1.PB.1.SE.allele_as' / 'haplo2.flnc.sorted.bam'}",
        "",
        "chr1",
        100,
        350,
        "200",
        "250",
        "#00EE00",
        str(report_dir / "PB.1.PB.1.SE.allele_as"),
    ]
    assert cwd == str(report_dir)
    assert runner is dummy_runner
    scripts = "\n".join(cmd[-1] for cmd, _ in dummy_runner.commands if isinstance(cmd, list) and cmd[:2] == ["bash", "-c"])
    assert "bed2gpe.pl" in scripts
    assert "-g 13 haplo1isosOut.bed" in scripts
    assert "genePredToGtf file haplo1isosOut.gpe haplo1isosOut.gtf" in scripts
    assert "samtools view -b -N haplo1.reads.txt" in scripts
    assert "samtools cat haplo1.flnc.sorted.bam haplo2.flnc.sorted.bam" in scripts
    assert "Generated allele-specific AS report plots" in capsys.readouterr().out


def test_report_samples_follow_the_strain_merged_layout(tmp_path):
    """BUG-160: after a default (strain-merged) run the outputs live in <project>/<strain>_all; the report must read
    them there, once per strain, and still fall back to per-row directories when only those exist."""
    from iflas.config import IflasConfig, Options, RefSection, Sample
    from iflas.report.plots import select_report_samples
    rows = [Sample(sample_name=n, project_name="demo", ref_strain="B73", strain="B73", data_processed_location=f"/r/{n}.fa")
            for n in ("s1", "s2")]
    cfg = IflasConfig(reference=[RefSection(ref_strain="B73", ref_genome="/g.fa", ref_gtf="/g.gtf")], options=Options(threads=2),
                      samples=rows, dirs={"out_dir": str(tmp_path / "out")})
    (tmp_path / "out" / "demo" / "s1").mkdir(parents=True)
    assert [s.sample_name for s in select_report_samples(cfg)] == ["s1", "s2"]          # only per-row dirs exist
    (tmp_path / "out" / "demo" / "B73_all").mkdir()
    assert [s.sample_name for s in select_report_samples(cfg)] == ["B73_all"]           # merged run: one report
    assert [s.sample_name for s in select_report_samples(cfg, ["s2"])] == ["B73_all"]
    assert [s.sample_name for s in select_report_samples(cfg, ["B73_all"])] == ["B73_all"]
    cfg.options.merge_data_from_same_strain = False
    assert [s.sample_name for s in select_report_samples(cfg)] == ["s1", "B73_all"]     # s2 has only the merged dir
    import pytest
    with pytest.raises(ValueError, match="unknown sample"):
        select_report_samples(cfg, ["nope"])
