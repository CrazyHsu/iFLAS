from iflas.stages.palen_as import (
    getASpairedIsoforms,
    getPalenAPA,
    getPalenAS,
    palen_as_sample,
)


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


def _write_pair_inputs(tmp_path):
    isoform_file = tmp_path / "hq.collapsed.bed12+"
    isoform_file.write_text(
        "\n".join(
            [
                _bed12plus("PB.1.1", [(100, 200), (300, 400)], "PB.1.1"),
                _bed12plus("PB.1.2", [(100, 200), (300, 400)], "PB.1.2"),
                _bed12plus("PB.1.3", [(100, 200), (500, 600)], "PB.1.3"),
                _bed12plus("PB.2.1", [(700, 900)], "PB.2.1", gene="geneY"),
            ]
        )
        + "\n"
    )
    group_file = tmp_path / "tofu.collapsed.group.txt"
    group_file.write_text(
        "\n".join(
            [
                "PB.1.1\tread1,read2",
                "PB.1.2\tread3",
                "PB.1.3\tread4,read5,read6",
                "PB.2.1\tread7,read8,read9",
            ]
        )
        + "\n"
    )
    return isoform_file, group_file


def _write_palen_base(tmp_path):
    base = tmp_path / "sample"
    hq_dir = base / "hqIsoforms"
    collapse_dir = base / "collapse"
    mapping_dir = base / "mapping"
    as_dir = base / "as_events" / "ordinary_as" / "LR"
    pa_dir = base / "as_events" / "pa"
    for d in (hq_dir, collapse_dir, mapping_dir, as_dir, pa_dir):
        d.mkdir(parents=True)

    reads_inc = [f"inc{i}" for i in range(30)]
    reads_exc = [f"exc{i}" for i in range(30)]
    (hq_dir / "hq.collapsed.bed12+").write_text(
        "\n".join(
            [
                _bed12plus("PB.1.1", [(100, 200), (300, 400)], "PB.1.1"),
                _bed12plus("PB.1.2", [(100, 200), (500, 600)], "PB.1.2"),
            ]
        )
        + "\n"
    )
    (collapse_dir / "tofu.collapsed.group.txt").write_text(
        "\n".join(
            [
                "PB.1.1\t" + ",".join(reads_inc),
                "PB.1.2\t" + ",".join(reads_exc),
            ]
        )
        + "\n"
    )
    (mapping_dir / "flnc.addCVandID.bed12+").write_text("")

    se_row = [
        "chr1",
        "200",
        "300",
        "PB.1:se",
        "0",
        "+",
        "200",
        "300",
        "0,0,0",
        "1",
        "100",
        "0",
        "geneX",
        "200",
        "300",
        "PB.1.1",
        "30",
        "PB.1.2",
        "30",
    ]
    (as_dir / "SE.confident.bed12+").write_text("\t".join(se_row) + "\n")
    for name in ("IR.confident.bed6+", "A3SS.confident.bed6+", "A5SS.confident.bed6+"):
        (as_dir / name).write_text("")

    (pa_dir / "paCluster.bed8+").write_text(
        "\n".join(
            [
                "\t".join(["chr1", "100", "120", ",".join(reads_inc), "30", "+", "0", "110", "PB.1"]),
                "\t".join(["chr1", "200", "220", ",".join(reads_exc), "30", "+", "0", "210", "PB.1"]),
            ]
        )
        + "\n"
    )
    palen = tmp_path / "polya_results.tsv"
    palen.write_text(
        "readname\tqc_tag\tpolya_length\n"
        + "\n".join([f"{read}\tPASS\t{10 + i}" for i, read in enumerate(reads_inc)])
        + "\n"
        + "\n".join([f"{read}\tPASS\t{200 + i}" for i, read in enumerate(reads_exc)])
        + "\nfailed\tFAIL\t999\n"
    )
    return base, palen


def test_get_as_paired_isoforms_filters_by_fl_count_and_merges_junctions(tmp_path):
    isoform_file, group_file = _write_pair_inputs(tmp_path)
    as_file = tmp_path / "A5SS.confident.bed6+"
    as_file.write_text(
        "\t".join(
            [
                "chr1",
                "200",
                "300",
                "PB.1:a5",
                "0",
                "+",
                "geneX",
                "PB.1.1,PB.1.2",
                "3",
                "PB.1.3",
                "3",
            ]
        )
        + "\n"
    )

    pairs = getASpairedIsoforms(
        str(as_file),
        str(group_file),
        str(isoform_file),
        asType="A5SS",
        filterByCount=2,
        mergeByJunc=True,
    )

    assert pairs == {"PB.1": [[["PB.1.1"], ["PB.1.3"], "PB.1:a5"]]}


def test_get_as_paired_isoforms_skips_blank_event_lines(tmp_path):
    isoform_file, group_file = _write_pair_inputs(tmp_path)
    as_file = tmp_path / "A5SS.confident.bed6+"
    as_file.write_text("\n\n")

    pairs = getASpairedIsoforms(
        str(as_file),
        str(group_file),
        str(isoform_file),
        asType="A5SS",
        filterByCount=0,
        mergeByJunc=True,
    )

    assert pairs == {}


def test_get_as_paired_isoforms_reads_se_inclusion_exclusion_columns(tmp_path):
    isoform_file, group_file = _write_pair_inputs(tmp_path)
    as_file = tmp_path / "SE.confident.bed12+"
    row = [
        "chr1",
        "200",
        "300",
        "PB.1:se",
        "0",
        "+",
        "200",
        "300",
        "0,0,0",
        "1",
        "100",
        "0",
        "geneX",
        "200",
        "300",
        "PB.1.1,PB.1.2",
        "3",
        "PB.1.3",
        "3",
    ]
    as_file.write_text("\t".join(row) + "\n")

    pairs = getASpairedIsoforms(
        str(as_file),
        str(group_file),
        str(isoform_file),
        asType="SE",
        filterByCount=0,
        mergeByJunc=False,
    )

    assert pairs == {"PB.1": [[["PB.1.1", "PB.1.2"], ["PB.1.3"], "PB.1:se"]]}


def test_get_palen_apa_reports_genes_with_significant_pa_length_shift(tmp_path):
    reads_a = [f"a{i}" for i in range(30)]
    reads_b = [f"b{i}" for i in range(30)]
    pa_cluster = tmp_path / "paCluster.bed8+"
    pa_cluster.write_text(
        "\n".join(
            [
                "\t".join(["chr1", "100", "120", ",".join(reads_a), "30", "+", "0", "110", "PB.1"]),
                "\t".join(["chr1", "200", "220", ",".join(reads_b), "30", "+", "0", "210", "PB.1"]),
            ]
        )
        + "\n"
    )
    flnc_reads_to_palen = {
        **{read: str(10 + i) for i, read in enumerate(reads_a)},
        **{read: str(200 + i) for i, read in enumerate(reads_b)},
    }
    out_file = tmp_path / "apaRelatedPalen.txt"

    getPalenAPA(str(pa_cluster), flnc_reads_to_palen, 5, str(out_file))

    fields = out_file.read_text().strip().split("\t")
    assert fields[:4] == ["PB.1", "+", "chr1_110;chr1_210", "30;30"]
    assert fields[4].split(";")[0].startswith("10.0,11.0,12.0")
    assert fields[4].split(";")[1].startswith("200.0,201.0,202.0")


def test_get_palen_as_writes_significant_pair_files(tmp_path, monkeypatch):
    reads_inc = [f"inc{i}" for i in range(30)]
    reads_exc = [f"exc{i}" for i in range(30)]
    isoform_file = tmp_path / "hq.collapsed.bed12+"
    isoform_file.write_text(
        "\n".join(
            [
                _bed12plus("PB.1.1", [(100, 200), (300, 400)], "PB.1.1"),
                _bed12plus("PB.1.2", [(100, 200), (500, 600)], "PB.1.2"),
            ]
        )
        + "\n"
    )
    reads_file = tmp_path / "flnc.addCVandID.bed12+"
    reads_file.write_text("")
    flnc_reads_to_palen = {
        **{read: str(10 + i) for i, read in enumerate(reads_inc)},
        **{read: str(200 + i) for i, read in enumerate(reads_exc)},
    }
    collapsed_trans_to_reads = {
        "PB.1.1": reads_inc,
        "PB.1.2": reads_exc,
    }
    as_pairs = {"SE": {"PB.1": [[["PB.1.1"], ["PB.1.2"], "PB.1:se"]]}}

    monkeypatch.chdir(tmp_path)
    getPalenAS(
        flnc_reads_to_palen,
        str(isoform_file),
        str(reads_file),
        collapsedTrans2reads=collapsed_trans_to_reads,
        asPairs=as_pairs,
        filterByCount=5,
        mergeByJunc=True,
    )

    sig_file = tmp_path / "mergeByJunc" / "SE.palenAndAS.sig.bed12+"
    sig_lines = sig_file.read_text().strip().splitlines()
    assert len(sig_lines) == 2
    assert sig_lines[0].split("\t")[12:] == [
        "geneX",
        "geneX",
        "PB.1.1",
        "PB.1.1+PB.1.2",
        "inclusion",
        ",".join(f"{10 + i}.0" for i in range(30)),
        "30",
    ]
    assert sig_lines[1].split("\t")[15:17] == ["PB.1.1+PB.1.2", "exclusion"]

    pair_file = tmp_path / "mergeByJunc" / "SE_PB.1_1.txt"
    pair_lines = pair_file.read_text().strip().splitlines()
    assert pair_lines[0] == "PB.1\tSE\tPB.1.1\tinc0\t10"
    assert pair_lines[-1] == "PB.1\tSE\tPB.1.2\texc29\t229"


def test_palen_as_sample_wires_legacy_file_contract(tmp_path):
    base, palen = _write_palen_base(tmp_path)

    out = palen_as_sample(str(base), [str(palen)], filter_by_count=5)

    sig_file = out["as_sig"]["SE"]
    sig_lines = sig_file.read_text().strip().splitlines()
    assert len(sig_lines) == 2
    assert sig_lines[0].split("\t")[15:17] == ["PB.1.1+PB.1.2", "inclusion"]
    assert out["as_sig"]["IR"].read_text() == ""

    apa_fields = out["apa"].read_text().strip().split("\t")
    assert apa_fields[:4] == ["PB.1", "+", "chr1_110;chr1_210", "30;30"]


class _NanopolishFakeRunner:
    def __init__(self):
        self.calls = []

    def run(self, cmd, **kwargs):
        import subprocess
        self.calls.append((cmd, kwargs))
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")


def test_polya_len_calling_follows_legacy_nanopolish_contract(tmp_path):
    from iflas.stages.palen_as import polya_len_calling

    fast5 = tmp_path / "fast5"
    fast5.mkdir()
    (fast5 / "a.fast5").write_text("")
    (fast5 / "sequencing_summary.txt").write_text("")
    runner = _NanopolishFakeRunner()

    out = polya_len_calling(str(fast5), "/p/rawFlnc.fq", "/m/flnc.mm2.sorted.bam", "/ref/genome.fa", tmp_path / "palenAS",
                            sequencing_summary=str(fast5 / "sequencing_summary.txt"), threads=4, runner=runner)

    assert out == str(tmp_path / "palenAS" / "polya_results.tsv")
    cmds = [c for c, _ in runner.calls]
    assert cmds[0] == f"nanopolish index --directory={fast5} --sequencing-summary={fast5 / 'sequencing_summary.txt'} /p/rawFlnc.fq"
    assert cmds[1] == "nanopolish polya --threads=4 --reads=/p/rawFlnc.fq --bam=/m/flnc.mm2.sorted.bam --genome=/ref/genome.fa > polya_results.tsv"
    assert "hist.R" in cmds[2] and "-p=polyaTailLength.pdf" in cmds[2] and runner.calls[2][1]["check"] is False
    assert all(kw["cwd"] == str(tmp_path / "palenAS") for _, kw in runner.calls)


def test_polya_len_calling_rejects_non_fast5_dir(tmp_path):
    import pytest
    from iflas.stages.palen_as import polya_len_calling

    (tmp_path / "reads.fq").write_text("")
    with pytest.raises(ValueError):
        polya_len_calling(str(tmp_path), "/p/rawFlnc.fq", "/m/bam", "/ref/genome.fa", tmp_path / "out", runner=_NanopolishFakeRunner())


def test_palen_as_sample_requires_palen_or_fast5(tmp_path):
    import pytest
    from iflas.stages.palen_as import palen_as_sample

    with pytest.raises(ValueError):
        palen_as_sample(str(tmp_path), [])
