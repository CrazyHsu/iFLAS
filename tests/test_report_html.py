import json
from pathlib import Path

from iflas.cli import main


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


def test_iflas_report_html_collects_existing_report_and_visual_as_outputs(tmp_path, capsys):
    out_dir = tmp_path / "out"
    cfg, samples = _write_report_config(tmp_path, out_dir)

    basic_dir = out_dir / "reports" / "demo_s1" / "basicStatistics"
    basic_dir.mkdir(parents=True)
    (basic_dir / "GC_of_raw_flnc_reads.png").write_text("fake png")

    visual_dir = out_dir / "demo" / "s1" / "visual_as" / "geneA"
    visual_dir.mkdir(parents=True)
    (visual_dir / "geneA.png").write_text("fake png")

    rc = main(["report", "--html", "-c", str(cfg), "-s", str(samples)])

    assert rc == 0
    html_dir = out_dir / "reports" / "html"
    main_page = html_dir / "iflas_report.html"
    basic_page = html_dir / "s1" / "basicStatistics.html"
    gene_page = html_dir / "s1" / "geneA.html"
    search_json = out_dir / "reports" / "json" / "search.json"

    assert main_page.is_file()
    assert basic_page.is_file()
    assert gene_page.is_file()
    assert "The GC content of flnc reads" in basic_page.read_text()
    assert "The isoform structure in geneA" in gene_page.read_text()
    assert (html_dir / "assets" / "generated" / "s1_basic_GC_of_raw_flnc_reads.png").is_file()
    assert (html_dir / "assets" / "generated" / "s1_gene_geneA_isoformStruc.png").is_file()

    data = json.loads(search_json.read_text())
    assert data["code"] == 0
    assert {item["title"] for item in data["data"]} == {
        "mainPage",
        "basic statistics page for s1",
        "single gene page of geneA for s1",
    }
    assert "Wrote HTML report" in capsys.readouterr().out
