from iflas.cli import main
from iflas.stages.asas import asas_sample, relationshipBetweenAlleleSpecifiAndAS


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


def _write_phased_partial(path):
    path.write_text(
        "# IsoPhase human-readable output\n"
        "haplotype\tPB.1.1\tPB.1.2\n"
        "hapREF\t20\t0\n"
        "hapALT\t0\t20\n"
    )


def _write_asas_base(tmp_path):
    base = tmp_path / "sample"
    hq_dir = base / "hqIsoforms"
    collapse_dir = base / "collapse"
    as_dir = base / "as_events" / "ordinary_as" / "LR"
    locus_dir = base / "asas" / "by_loci" / "PB.1_size20"
    for directory in (hq_dir, collapse_dir, as_dir, locus_dir):
        directory.mkdir(parents=True)

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
                "PB.1.1\t" + ",".join(f"inc{i}" for i in range(20)),
                "PB.1.2\t" + ",".join(f"exc{i}" for i in range(20)),
            ]
        )
        + "\n"
    )
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
        "20",
        "PB.1.2",
        "20",
    ]
    (as_dir / "SE.confident.bed12+").write_text("\t".join(se_row) + "\n")
    for name in ("IR.confident.bed6+", "A5SS.confident.bed6+", "A3SS.confident.bed6+"):
        (as_dir / name).write_text("")
    _write_phased_partial(locus_dir / "phased.partial.cleaned.human_readable.txt")
    return base


def test_relationship_between_allele_specific_and_as_detects_partial_signal(tmp_path, monkeypatch):
    _write_phased_partial(tmp_path / "phased.partial.cleaned.human_readable.txt")
    monkeypatch.chdir(tmp_path)

    result = relationshipBetweenAlleleSpecifiAndAS(
        partial=True,
        nopartial=False,
        isoPairs={"PB.1:se": [["PB.1.1"], ["PB.1.2"]]},
    )

    assert result == {
        "partial": {
            "PB.1:se": {
                "hapREF": "PB.1.1",
                "hapALT": "PB.1.2",
            }
        },
        "nopartial": {},
    }


def test_asas_sample_wires_existing_phased_loci_to_legacy_output(tmp_path):
    base = _write_asas_base(tmp_path)

    out = asas_sample(str(base))

    lines = out["partial"].read_text().strip().splitlines()
    assert lines == [
        "\t".join(["PB.1", "SE", "PB.1", "PB.1:se", "hapREF", "PB.1.1", "hapALT", "PB.1.2"])
    ]


def test_cli_asas_runs_file_level_stage(tmp_path):
    base = _write_asas_base(tmp_path)

    rc = main(["asas", "-b", str(base)])

    assert rc == 0
    assert (base / "asas" / "partialAsRelatedHaplotype.txt").read_text().strip().split("\t")[:4] == [
        "PB.1",
        "SE",
        "PB.1",
        "PB.1:se",
    ]


class _PhaserFakeRunner:
    def __init__(self):
        self.calls = []

    def run(self, cmd, **kwargs):
        self.calls.append((cmd, kwargs))
        import subprocess
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")


def test_make_abundance_file_is_cupcake_read_stat(tmp_path):
    from iflas.stages.asas import makeAbundanceFile

    group = tmp_path / "tofu.collapsed.group.txt"
    group.write_text("PB.1.1\tr1,r2\nPB.1.2\tr3\n")
    makeAbundanceFile(str(group), str(tmp_path / "read_stat.txt"))
    assert (tmp_path / "read_stat.txt").read_text() == (
        "id\tlength\tis_fl\tstat\tpbid\nr1\tNA\tY\tunique\tPB.1.1\nr2\tNA\tY\tunique\tPB.1.1\nr3\tNA\tY\tunique\tPB.1.2\n")


def test_select_loci_to_phase_routes_to_cupcake_env_and_detects_fastq(tmp_path):
    from iflas.stages.asas import select_loci_to_phase

    flnc = tmp_path / "flnc.fq"
    flnc.write_text("@r1\nACGT\n+\nIIII\n")
    runner = _PhaserFakeRunner()
    loci = select_loci_to_phase("/ref/genome.fa", str(flnc), "/c/tofu.collapsed.good.gff", "/a/read_stat.txt", tmp_path,
                                min_coverage=10, runner=runner)
    cmd, kwargs = runner.calls[0]
    assert cmd == f"python $(which select_loci_to_phase.py) /ref/genome.fa {flnc} /c/tofu.collapsed.good.gff /a/read_stat.txt -c 10 --fq"
    assert kwargs["conda_env"] == "iflas-cupcake" and kwargs["cwd"] == str(tmp_path)
    assert loci == []


def test_run_phaser_follows_legacy_command_contract(tmp_path):
    from iflas.stages.asas import run_phaser

    locus = tmp_path / "by_loci" / "PB.3_size12"
    locus.mkdir(parents=True)
    (locus / "config").write_text("pbid=PB.3\nref_chr=chr2\nref_strand=-\nref_start=100\nref_end=900\n")
    runner = _PhaserFakeRunner()

    assert run_phaser(locus, runner=runner) is True

    cmds = [c for c, _ in runner.calls]
    assert cmds[0] == "minimap2 -ax splice fake.fasta ccs.fastq > ccs.sam 2> /dev/null"
    assert cmds[1].startswith("samtools sort ccs.sam > ccs.sorted.bam") and "samtools mpileup --min-BQ 13 -f fake.fasta -s ccs.sorted.bam" in cmds[1]
    assert cmds[2] == ("run_phaser.py ccs.fastq ccs.sam ccs.mpileup fake.read_stat.txt fake.mapping.txt --partial_ok --strand - "
                       "-o phased.partial -n 2 --bhFDR 0.01 -p 0.01 > partial.stdout.txt 2>&1")
    assert runner.calls[2][1]["conda_env"] == "iflas-cupcake" and runner.calls[2][1]["check"] is False
    assert all(kw["cwd"] == str(locus) for _, kw in runner.calls)


_VCF = """##fileformat=VCFv4.2
##contig=<ID=chr2,length=100000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr2\t40100\t.\tC\tT\t50\tPASS\tDP=30
chr2\t40300\t.\tA\tC\t50\tPASS\tDP=30
chr2\t40500\t.\tA\tT\t50\tPASS\tDP=30
"""


def test_vcf_haplo_and_assign_haplo_label_isophase_haplotypes(tmp_path):
    from iflas.stages.asas import assignHaplo, vcfHaplo

    vcf = tmp_path / "phased.partial.cleaned.vcf"
    vcf.write_text(_VCF)
    labels = vcfHaplo(str(vcf))
    assert labels == {"CAA": "REF", "TCT": "ALT"}
    assert assignHaplo(["CAA", "TCT", "CAT"], labels) == {"CAA": "REF", "TCT": "ALT", "CAT": "REF"}   # closest match


def test_generate_regions_tiles_the_fasta_index(tmp_path):
    from iflas.stages.asas import generate_regions

    (tmp_path / "g.fa.fai").write_text("chr1\t2500\t0\t60\t61\nchr2\t900\t0\t60\t61\n")
    generate_regions(str(tmp_path / "g.fa"), 1000, str(tmp_path / "regions.txt"))
    assert (tmp_path / "regions.txt").read_text().splitlines() == ["chr1:0-1000", "chr1:1000-2000", "chr1:2000-2500", "chr2:0-900"]


def test_merge_ase_from_salmon_out_joins_fl_and_ngs_support(tmp_path):
    import pandas as pd
    from iflas.stages.asas import mergeASEfromSalmonOut

    ref = tmp_path / "ref.quant.sf"
    alt = tmp_path / "alt.quant.sf"
    pd.DataFrame({"Name": ["PB.48_size58_REF"], "Length": [700], "EffectiveLength": [500.0], "TPM": [120.5], "NumReads": [40.0]}).to_csv(ref, sep="\t", index=False)
    pd.DataFrame({"Name": ["PB.48_size58_ALT"], "Length": [700], "EffectiveLength": [500.0], "TPM": [30.25], "NumReads": [10.0]}).to_csv(alt, sep="\t", index=False)
    haplo2count = {"PB.48_size58_REF": [21, "asasA"], "PB.48_size58_ALT": [5, "asasA"]}

    mergeASEfromSalmonOut(str(ref), str(alt), haplo2count, str(tmp_path / "merged.txt"))

    assert (tmp_path / "merged.txt").read_text() == "PB.48\tasasA\t21.0\t5.0\t120.5\t40.0\t30.25\t10.0\n"


def test_bug180_haplotype_without_reads_of_either_side_does_not_crash(tmp_path, monkeypatch):
    """IsoPhase tables can hold more than two haplotypes (partial ones); one without a read of either side of an event
    gave chi2_contingency a zero row ("expected frequencies has a zero element") and aborted the whole asas stage."""
    (tmp_path / "phased.partial.cleaned.human_readable.txt").write_text(
        "# IsoPhase\nhaplotype\tPB.1.1\tPB.1.2\tPB.1.3\nhapA\t30\t0\t2\nhapB\t0\t30\t0\nhapC\t0\t0\t9\n")
    monkeypatch.chdir(tmp_path)
    result = relationshipBetweenAlleleSpecifiAndAS(partial=True, isoPairs={"PB.1:se": [["PB.1.1"], ["PB.1.2"]]})
    assert result["partial"] == {"PB.1:se": {"hapA": "PB.1.1", "hapB": "PB.1.2"}}


def test_bug201_phased_loci_need_an_snp_the_ngs_reads_confirm(tmp_path):
    """BUG-201 (adopted 2026-10-16): IsoPhase reads isoform-specific alignment differences at splice sites as
    SNPs, so the association test called AS events allele-specific without any variant; a locus is kept only when
    freebayes also calls one of its SNPs on the NGS reads, and ase uses the same loci."""
    from iflas.stages.asas import ngs_confirmed_loci

    header = "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n##INFO=<ID=TYPE,Number=A,Type=String,Description=\"type\">\n" \
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    loci = []
    for name, pos in (("PB.1_size20", 100), ("PB.2_size20", 500)):
        locus = tmp_path / "by_loci" / name
        locus.mkdir(parents=True)
        (locus / "phased.partial.cleaned.human_readable.txt").write_text("Associated VCF file\nhaplotype\tPB.1.1\nA\t5\n")
        (locus / "phased.partial.cleaned.vcf").write_text(header + f"chr1\t{pos}\t.\tA\tG\t.\t.\tTYPE=snp\n")
        loci.append(locus)

    class FreebayesConfirmsTheFirst:
        def run(self, cmd, **kwargs):
            (tmp_path / "phased_snps.ngs.vcf").write_text(header + "chr1\t100\t.\tA\tG\t50\t.\tTYPE=snp\n")

    bam = tmp_path / "tmp.bam"
    bam.write_bytes(b"")
    kept = ngs_confirmed_loci(loci, "genome.fa", str(bam), tmp_path, FreebayesConfirmsTheFirst())
    assert [locus.name for locus in kept] == ["PB.1_size20"]
    assert (tmp_path / "ngs_confirmed_loci.txt").read_text() == "PB.1_size20\n"
    assert [locus.name for locus in ngs_confirmed_loci(loci, "genome.fa", str(tmp_path / "none.bam"), tmp_path, None)] == \
        ["PB.1_size20", "PB.2_size20"]                                       # no NGS alignments: every locus kept
