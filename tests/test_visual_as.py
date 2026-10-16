"""Tests for the migrated file-level visual_as pyGenomeTracks bridge."""

from types import SimpleNamespace

from iflas.cli import main
from iflas.stages.visual_as import build_visual_as_tracks, visual_as_merged


def _bed12plus(name, exons, group="PB.1.1", gene="geneA"):
    start = exons[0][0]
    end = exons[-1][1]
    sizes = ",".join(str(e - s) for s, e in exons)
    rel_starts = ",".join(str(s - start) for s, e in exons)
    fields = [
        "chr1",
        start,
        end,
        name,
        0,
        "+",
        start,
        end,
        "0,0,0",
        len(exons),
        sizes,
        rel_starts,
        group.rsplit(".", 1)[0],     # col 13: cupcake gene id (PB.1)
        gene,                        # col 14: gene name
        group,                       # col 15: group key
    ]
    return "\t".join(map(str, fields))


def _gpe(transcript, exons, gene="geneA"):
    starts = ",".join(str(s) for s, _ in exons) + ","
    ends = ",".join(str(e) for _, e in exons) + ","
    return "\t".join(
        map(
            str,
            [
                transcript,
                "chr1",
                "+",
                exons[0][0],
                exons[-1][1],
                exons[0][0],
                exons[-1][1],
                len(exons),
                starts,
                ends,
                0,
                gene,
                "cmpl",
                "cmpl",
                ",".join(["0"] * len(exons)) + ",",
            ],
        )
    )


def test_build_visual_as_tracks_writes_beds_and_pygenometracks_config(tmp_path):
    base = tmp_path / "sample"
    hq_dir = base / "hqIsoforms"
    hq_dir.mkdir(parents=True)
    isoform_bed = hq_dir / "hq.collapsed.bed12+"
    isoform_bed.write_text(
        _bed12plus("PB.1.1", [(100, 180), (260, 320)]) + "\n"
        + _bed12plus("PB.1.2", [(100, 320)]) + "\n"
    )
    ref_gpe = tmp_path / "ref.gpe"
    ref_gpe.write_text(_gpe("geneA.1", [(90, 180), (260, 330)]) + "\n")

    out = build_visual_as_tracks(base, ref_gpe, genes=["geneA"])

    assert out["region"] == "chr1:90-330"
    assert out["gene_dir"].endswith("visual_as/geneA")
    assert "file = " in (tmp_path / "sample" / "visual_as" / "geneA" / "tracks.ini").read_text()
    assert (tmp_path / "sample" / "visual_as" / "geneA" / "geneA.reference.bed").read_text().startswith(
        "chr1\t90\t330\tgeneA.1"
    )
    iso_bed = (tmp_path / "sample" / "visual_as" / "geneA" / "geneA.iflas_isoforms.bed").read_text()
    assert iso_bed.count("\n") == 2 and all(len(l.split("\t")) == 12 for l in iso_bed.splitlines())   # plain bed12 for pyGenomeTracks


def test_iflas_visual_as_cli_runs_pygenometracks_contract(tmp_path, monkeypatch, capsys):
    base = tmp_path / "sample"
    hq_dir = base / "hqIsoforms"
    hq_dir.mkdir(parents=True)
    isoform_bed = hq_dir / "hq.collapsed.bed12+"
    isoform_bed.write_text(_bed12plus("PB.1.1", [(100, 180), (260, 320)]) + "\n")
    ref_gpe = tmp_path / "ref.gpe"
    ref_gpe.write_text(_gpe("geneA.1", [(90, 180), (260, 330)]) + "\n")

    class FakeRunner:
        def __init__(self):
            self.calls = []

        def run(self, cmd, **kwargs):
            self.calls.append((cmd, kwargs))
            return SimpleNamespace(returncode=0, stdout="", stderr="")

    runner = FakeRunner()

    import iflas.stages.visual_as

    monkeypatch.setattr(iflas.stages.visual_as, "ToolRunner", lambda: runner)
    rendered = []                                   # G4: the CLI renders in-process through pyGenomeTracks' entry point
    monkeypatch.setattr(iflas.stages.visual_as, "render_all", lambda jobs, threads: rendered.extend(jobs))
    rc = main(["visual_as", "-b", str(base), "-r", str(ref_gpe), "--gene", "geneA"])

    assert rc == 0
    out = capsys.readouterr().out
    assert "Wrote visual_as outputs" in out
    assert runner.calls == []
    assert [(str(t), r, str(p)) for t, r, p in rendered] == [(str(base / "visual_as" / "geneA" / "tracks.ini"), "chr1:90-330", str(base / "visual_as" / "geneA" / "geneA.pdf"))]


def test_iflas_visual_as_cli_gene_file_skips_comments_and_blank_lines(tmp_path, monkeypatch, capsys):
    base = tmp_path / "sample"
    hq_dir = base / "hqIsoforms"
    hq_dir.mkdir(parents=True)
    isoform_bed = hq_dir / "hq.collapsed.bed12+"
    isoform_bed.write_text(
        _bed12plus("PB.1.1", [(100, 180), (260, 320)], gene="geneA") + "\n"
        + _bed12plus("PB.2.1", [(500, 560), (640, 700)], group="PB.2.1", gene="geneB") + "\n"
    )
    ref_gpe = tmp_path / "ref.gpe"
    ref_gpe.write_text(
        _gpe("geneA.1", [(90, 180), (260, 330)], gene="geneA") + "\n"
        + _gpe("geneB.1", [(490, 560), (640, 710)], gene="geneB") + "\n"
    )
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("\n# target genes\n geneB\tpreferred\n")

    class FakeRunner:
        def __init__(self):
            self.calls = []

        def run(self, cmd, **kwargs):
            self.calls.append((cmd, kwargs))
            return SimpleNamespace(returncode=0, stdout="", stderr="")

    runner = FakeRunner()

    import iflas.stages.visual_as

    monkeypatch.setattr(iflas.stages.visual_as, "ToolRunner", lambda: runner)
    rendered = []                                   # G4: the CLI renders in-process through pyGenomeTracks' entry point
    monkeypatch.setattr(iflas.stages.visual_as, "render_all", lambda jobs, threads: rendered.extend(jobs))
    rc = main(["visual_as", "-b", str(base), "-r", str(ref_gpe), "--gene-file", str(gene_file)])

    assert rc == 0
    out = capsys.readouterr().out
    assert "geneB" in out
    assert runner.calls == []
    assert [(str(t), r, str(p)) for t, r, p in rendered] == [(str(base / "visual_as" / "geneB" / "tracks.ini"), "chr1:490-710", str(base / "visual_as" / "geneB" / "geneB.pdf"))]


class _NgsFakeRunner:
    """Answers `samtools depth` with a fixed 5-position profile and records every command."""

    def __init__(self):
        self.calls = []

    def run(self, cmd, **kwargs):
        self.calls.append((cmd, kwargs))
        if cmd[:2] == ["samtools", "depth"]:
            return SimpleNamespace(returncode=0, stdout="chr1\t91\t3\nchr1\t92\t3\nchr1\t93\t5\nchr1\t94\t0\nchr1\t95\t0\n", stderr="")
        return SimpleNamespace(returncode=0, stdout="", stderr="")


def _junctions_bed(path):
    # regtools bed12: intron = chromStart + size0 .. chromStart + blockStart1
    path.write_text("chr1\t170\t270\tJUNC1\t15\t+\t170\t270\t255,0,0\t2\t10,10\t0,90\n"      # 180-260 inside geneA
                    "chr1\t900\t1000\tJUNC2\t20\t+\t900\t1000\t255,0,0\t2\t10,10\t0,90\n")  # elsewhere


def test_build_visual_as_tracks_adds_ngs_coverage_and_junction_tracks(tmp_path):
    base = tmp_path / "sampleA"
    (base / "hqIsoforms").mkdir(parents=True)
    (base / "hqIsoforms" / "hq.collapsed.bed12+").write_text(_bed12plus("PB.1.1", [(100, 180), (260, 320)]) + "\n")
    ref_gpe = tmp_path / "ref.gpe"
    ref_gpe.write_text(_gpe("geneA.1", [(90, 180), (260, 330)]) + "\n")
    junctions = tmp_path / "junctions.bed"
    _junctions_bed(junctions)
    runner = _NgsFakeRunner()

    out = build_visual_as_tracks(base, ref_gpe, genes=["geneA"], runner=runner, run_plot=True,
                                 ngs_bams=[tmp_path / "repeat0.sorted.bam"], ngs_junctions=junctions)

    gene_dir = base / "visual_as" / "geneA"
    assert runner.calls[0][0] == ["samtools", "depth", "-a", "-r", "chr1:91-330", str(tmp_path / "repeat0.sorted.bam")]
    # runs of equal depth become one interval (per-base lines made long genes' tracks huge and slow to draw)
    assert (gene_dir / "repeat0.sampleA.bedgraph").read_text() == "chr1\t90\t92\t3\nchr1\t92\t93\t5\nchr1\t93\t95\t0\n"
    assert (gene_dir / "sampleA.junctions.links").read_text() == "chr1\t179\t180\tchr1\t260\t261\t15\n"   # only the in-region junction
    ini = (gene_dir / "tracks.ini").read_text()
    assert "[reference]" in ini and "[iflas_isoforms]" in ini and "[coverage_repeat0_sampleA]" in ini and "[junctions_sampleA]" in ini
    assert out["ngs_tracks"] == 2 and out["region"] == "chr1:90-330"
    assert runner.calls[-1][0][:3] == ["pyGenomeTracks", "--tracks", str(gene_dir / "tracks.ini")]


def test_build_visual_as_tracks_novel_cluster_without_reference_model(tmp_path):
    base = tmp_path / "sampleA"
    (base / "hqIsoforms").mkdir(parents=True)
    (base / "hqIsoforms" / "hq.collapsed.bed12+").write_text(_bed12plus("PB.9.1", [(5000, 5100), (5300, 5400)], gene="chr1:+:3") + "\n")
    ref_gpe = tmp_path / "ref.gpe"
    ref_gpe.write_text(_gpe("geneA.1", [(90, 180), (260, 330)]) + "\n")

    out = build_visual_as_tracks(base, ref_gpe, genes=["chr1:+:3"])

    assert out["region"] == "chr1:5000-5400" and out["reference_bed"] == ""
    assert "[reference]" not in (base / "visual_as" / "chr1:+:3" / "tracks.ini").read_text()


def test_visual_as_merged_stacks_samples_under_one_reference(tmp_path):
    ref_gpe = tmp_path / "ref.gpe"
    ref_gpe.write_text(_gpe("geneA.1", [(90, 180), (260, 330)]) + "\n")
    bases = []
    for name, exons in (("s1", [(100, 180), (260, 320)]), ("s2", [(80, 180), (260, 340)])):
        base = tmp_path / "out" / "demo" / name
        (base / "hqIsoforms").mkdir(parents=True)
        (base / "hqIsoforms" / "hq.collapsed.bed12+").write_text(_bed12plus("PB.1.1", exons) + "\n")
        bases.append(base)
    runner = _NgsFakeRunner()

    out = visual_as_merged(bases, ref_gpe, genes="geneA", runner=runner, run_plot=True)

    gene_dir = tmp_path / "out" / "isoViewer_sample_merged" / "geneA"
    ini = (gene_dir / "tracks.ini").read_text()
    assert "[isoforms_s1]" in ini and "[isoforms_s2]" in ini and ini.count("file_type = bed\n") == 3
    assert out == [{"gene": "geneA", "gene_dir": str(gene_dir), "tracks": str(gene_dir / "tracks.ini"),
                    "pdf": str(gene_dir / "geneA.pdf"), "region": "chr1:80-340", "samples": ["s1", "s2"]}]
    assert runner.calls[-1][0] == ["pyGenomeTracks", "--tracks", str(gene_dir / "tracks.ini"), "--region", "chr1:80-340",
                                   "--outFileName", str(gene_dir / "geneA.pdf")]


def test_junction_index_selects_the_same_links_as_a_full_scan(tmp_path):
    """G4: genes pick their NGS junctions from a per-sample index; the links written must equal a full file scan
    (same junctions, same order) for regions anywhere along the chromosome."""
    import random
    from iflas.stages.visual_as import junction_index, write_ngs_tracks
    rng = random.Random(1)
    lines = []
    for k in range(300):
        cs = rng.randint(0, 50000)
        intron = rng.randint(60, 3000)
        chrom = rng.choice(["chr1", "chr2"])
        lines.append("\t".join(map(str, [chrom, cs, cs + intron + 20, f"j{k}", rng.randint(1, 50), "+", cs, cs + intron + 20,
                                         "0,0,0", 2, "10,10", f"0,{intron + 10}"])))
    bed = tmp_path / "junctions.bed"
    bed.write_text("\n".join(lines) + "\n")
    index = junction_index(str(bed))
    for k, (chrom, start, end) in enumerate([("chr1", 0, 60000), ("chr1", 20000, 21000), ("chr2", 49000, 52000), ("chr3", 0, 10)]):
        expected = []
        for line in lines:
            c = line.split("\t")
            sizes, starts = c[10].split(","), c[11].split(",")
            s, e = int(c[1]) + int(sizes[0]), int(c[1]) + int(starts[1])
            if c[0] == chrom and e > start and s < end:
                expected.append(f"{chrom}\t{s - 1}\t{s}\t{chrom}\t{e}\t{e + 1}\t{c[4]}\n")
        gene_dir = tmp_path / f"g{k}"
        gene_dir.mkdir()
        write_ngs_tracks(gene_dir, "S", chrom, start, end, [], index, runner=None)
        assert (gene_dir / "S.junctions.links").read_text() == "".join(expected)
