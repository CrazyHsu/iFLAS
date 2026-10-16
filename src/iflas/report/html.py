"""HTML report aggregation for already-generated iFLAS outputs.

This is the file-level, no-external-tool slice of legacy ``generateHtml.py``.
It discovers existing plot outputs, materializes them under ``reports/html``,
and writes simple HTML pages plus the legacy ``reports/json/search.json`` index.
Plot generation itself still lives in the individual migrated stages or the
legacy Rscript-backed helpers.
"""

from __future__ import annotations

import html
import json
import os
import re
import shutil
from pathlib import Path
from typing import Iterable

from ..config import IflasConfig, Sample

_FIGURE_SUFFIXES = (".png", ".jpg", ".jpeg", ".gif", ".svg", ".pdf")

_BASIC_FIGURES = (
    ("readsCorrectResult", "The evaluation of reads before and after correction"),
    ("GC_of_raw_flnc_reads", "The GC content of flnc reads"),
    ("GC_across_raw_flnc_read", "The GC content across flnc reads"),
    ("LengthDistribution.curve", "The length distribution of full-length reads"),
    ("LengthDistribution.box", "The length distribution of full-length reads"),
    ("supportedByRNAseq", "The junctions of full-length reads supported by NGS reads"),
)


def _slug(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return value.strip("._") or "item"


def _find_figure(directory: Path, stem: str) -> Path | None:
    for suffix in _FIGURE_SUFFIXES:
        candidate = directory / f"{stem}{suffix}"
        if candidate.is_file():
            return candidate
    return None


def _convert_pdf_to_png(source: Path, dest: Path, page_index: int = 0) -> None:
    import fitz

    with fitz.open(source) as doc:
        page = doc.load_page(page_index)
        pix = page.get_pixmap(matrix=fitz.Matrix(2.0, 2.0))
        pix.save(dest)


def _materialize_figure(source: Path, html_dir: Path, name: str) -> Path:
    generated_dir = html_dir / "assets" / "generated"
    generated_dir.mkdir(parents=True, exist_ok=True)
    suffix = ".png" if source.suffix.lower() == ".pdf" else source.suffix.lower()
    dest = generated_dir / f"{_slug(name)}{suffix}"
    if source.suffix.lower() == ".pdf":
        _convert_pdf_to_png(source, dest)
    else:
        shutil.copy2(source, dest)
    return dest


def _rel(path: Path, from_dir: Path) -> str:
    return os.path.relpath(Path(path).resolve(), from_dir.resolve()).replace(os.sep, "/")


def _write_html(path: Path, title: str, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "<!doctype html>",
                '<html lang="en">',
                "<head>",
                '  <meta charset="utf-8">',
                '  <meta name="viewport" content="width=device-width, initial-scale=1">',
                f"  <title>{html.escape(title)}</title>",
                "  <style>",
                "    body { font-family: Arial, sans-serif; margin: 2rem; line-height: 1.5; }",
                "    nav a { margin-right: 1rem; }",
                "    img { max-width: 100%; height: auto; border: 1px solid #ddd; }",
                "    section { margin: 1.5rem 0; }",
                "  </style>",
                "</head>",
                "<body>",
                body,
                "</body>",
                "</html>",
                "",
            ]
        )
    )


def _collect_basic_statistics(report_dir: Path, project_name: str, sample_name: str) -> list[dict]:
    basic_dir = report_dir / "basicStatistics"
    as_pattern_dir = report_dir / "asPattern"
    figures: list[dict] = []
    if basic_dir.is_dir():
        for stem, title in _BASIC_FIGURES:
            source = _find_figure(basic_dir, stem)
            if source:
                figures.append({"key": stem, "title": title, "source": source})
    if as_pattern_dir.is_dir():
        for stem, title in (
            (f"{project_name}_{sample_name}.AS_annotation", "The alternative splicing(AS) summary"),
            (f"{project_name}_{sample_name}.AS_spliceSite", "The splice-site of alternative splicing(AS)"),
        ):
            source = _find_figure(as_pattern_dir, stem)
            if source:
                figures.append({"key": stem, "title": title, "source": source})
    return figures


def _collect_gene_figures(base_dir: Path) -> dict[str, dict]:
    genes: dict[str, dict] = {}
    for root_name in ("visual_as", "isoViewer"):
        root = base_dir / root_name
        if not root.is_dir():
            continue
        for gene_dir in sorted(path for path in root.iterdir() if path.is_dir()):
            gene = gene_dir.name
            source = _find_figure(gene_dir, gene)
            if not source:
                matches = [p for p in sorted(gene_dir.iterdir()) if p.suffix.lower() in _FIGURE_SUFFIXES]
                source = matches[0] if matches else None
            if source and gene not in genes:
                genes[gene] = {"isoformStruc": {"title": f"The isoform structure in {gene}", "source": source}}
    return genes


def collect_report_results(config: IflasConfig, samples: Iterable[Sample] | None = None) -> dict:
    """Discover report-ready files using the legacy output layout."""
    selected = list(samples if samples is not None else config.samples)
    out_dir = Path(config.dirs.out_dir).expanduser().resolve()
    result: dict = {}
    for sample in selected:
        project_name = sample.project_name
        sample_name = sample.sample_name
        base_dir = out_dir / project_name / sample_name
        report_dir = out_dir / "reports" / f"{project_name}_{sample_name}"
        result[sample_name] = {
            "basicStatistics": _collect_basic_statistics(report_dir, project_name, sample_name),
            "genes": _collect_gene_figures(base_dir),
        }
    return result


def _prepare_assets(report: dict, html_dir: Path) -> None:
    for sample_name, sample_report in report.items():
        for figure in sample_report.get("basicStatistics", []):
            figure["asset"] = _materialize_figure(
                figure["source"], html_dir, f"{sample_name}_basic_{figure['key']}"
            )
        for gene, gene_report in sample_report.get("genes", {}).items():
            for key, figure in gene_report.items():
                figure["asset"] = _materialize_figure(
                    figure["source"], html_dir, f"{sample_name}_gene_{gene}_{key}"
                )


def _main_page(report: dict, html_dir: Path) -> Path:
    rows = ["<h1>iFLAS report</h1>", "<p>Results of each sample.</p>", "<ul>"]
    for sample_name, sample_report in report.items():
        rows.append(f"  <li>{html.escape(sample_name)}<ul>")
        if sample_report.get("basicStatistics"):
            rows.append(
                f'    <li><a href="{html.escape(sample_name)}/basicStatistics.html">Basic Statistics</a></li>'
            )
        for gene in sample_report.get("genes", {}):
            rows.append(
                f'    <li><a href="{html.escape(sample_name)}/{html.escape(gene)}.html">'
                f"{html.escape(gene)}</a></li>"
            )
        rows.append("  </ul></li>")
    rows.append("</ul>")
    path = html_dir / "iflas_report.html"
    _write_html(path, "iFLAS report", "\n".join(rows))
    return path


def _sample_basic_page(sample_name: str, sample_report: dict, html_dir: Path) -> Path | None:
    figures = sample_report.get("basicStatistics", [])
    if not figures:
        return None
    page_dir = html_dir / sample_name
    rows = [
        '<nav><a href="../iflas_report.html">iFLAS</a></nav>',
        f"<h1>Basic statistics in sample {html.escape(sample_name)}</h1>",
    ]
    for figure in figures:
        src = _rel(figure["asset"], page_dir)
        title = html.escape(figure["title"])
        rows.extend(
            [
                "<section>",
                f"  <h3>{title}</h3>",
                f'  <img src="{html.escape(src)}" alt="{title}">',
                "</section>",
            ]
        )
    path = page_dir / "basicStatistics.html"
    _write_html(path, f"Basic statistics in sample {sample_name}", "\n".join(rows))
    return path


def _sample_gene_pages(sample_name: str, sample_report: dict, html_dir: Path) -> list[Path]:
    paths = []
    page_dir = html_dir / sample_name
    for gene, gene_report in sample_report.get("genes", {}).items():
        rows = [
            '<nav><a href="../iflas_report.html">iFLAS</a></nav>',
            f"<h1>Detailed information in gene {html.escape(gene)}</h1>",
        ]
        for figure in gene_report.values():
            src = _rel(figure["asset"], page_dir)
            title = html.escape(figure["title"])
            rows.extend(
                [
                    "<section>",
                    f"  <h2>{title}</h2>",
                    f'  <img src="{html.escape(src)}" alt="{title}">',
                    "</section>",
                ]
            )
        path = page_dir / f"{gene}.html"
        _write_html(path, f"Detailed information in gene {gene}", "\n".join(rows))
        paths.append(path)
    return paths


def generate_html_report(
    config: IflasConfig,
    *,
    sample_names: Iterable[str] | None = None,
) -> dict:
    """Generate ``reports/html`` and ``reports/json/search.json`` from existing outputs."""
    from .plots import select_report_samples

    selected = select_report_samples(config, sample_names)

    out_dir = Path(config.dirs.out_dir).expanduser().resolve()
    html_dir = out_dir / "reports" / "html"
    html_dir.mkdir(parents=True, exist_ok=True)
    report = collect_report_results(config, selected)
    _prepare_assets(report, html_dir)

    data = []
    main = _main_page(report, html_dir)
    data.append({"title": "mainPage", "url": str(main)})
    for sample_name, sample_report in report.items():
        basic = _sample_basic_page(sample_name, sample_report, html_dir)
        if basic:
            data.append({"title": f"basic statistics page for {sample_name}", "url": str(basic)})
        for gene_page in _sample_gene_pages(sample_name, sample_report, html_dir):
            gene = gene_page.stem
            data.append({"title": f"single gene page of {gene} for {sample_name}", "url": str(gene_page)})

    json_dir = out_dir / "reports" / "json"
    json_dir.mkdir(parents=True, exist_ok=True)
    json_file = json_dir / "search.json"
    json_file.write_text(json.dumps({"code": 0, "data": data}, indent=4))
    return {"html_dir": str(html_dir), "main_page": str(main), "json": str(json_file), "entries": data}
