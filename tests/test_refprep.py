"""End-to-end test for iflas.refprep.init_ref using a tiny synthetic genome.

Exercises real tools from the env (samtools faidx, gtfToGenePred, minimap2 -d).
"""

import os
import re
import shutil

import pytest

from iflas.cli import main
from iflas.refprep import ids_with_colons, init_ref

V2 = os.path.join(os.path.dirname(__file__), "..", "test_data", "synthetic_v2")


def _make_genome(tmp_path):
    fa = tmp_path / "genome.fa"
    fa.write_text(">chr1\n" + "ACGT" * 50 + "\n>chr2\n" + "TGCA" * 25 + "\n")  # 200bp, 100bp
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr1\tiflas\texon\t1\t50\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\tiflas\texon\t101\t150\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    )
    return str(fa), str(gtf)


def test_init_ref_end_to_end(tmp_path):
    fa, gtf = _make_genome(tmp_path)
    res = init_ref(fa, gtf=gtf, outdir=str(tmp_path), build_minimap2=True, build_hisat2=False)

    # chrom sizes correct (the legacy derivation was broken; BUG-110)
    sizes = dict(line.split("\t") for line in open(res["ref_size"]).read().splitlines())
    assert sizes["chr1"] == "200" and sizes["chr2"] == "100"

    # gpe + bed produced
    assert os.path.getsize(res["ref_gpe"]) > 0
    bed_lines = open(res["ref_bed"]).read().strip().splitlines()
    assert bed_lines and bed_lines[0].split("\t")[0] == "chr1"
    assert bed_lines[0].split("\t")[-1] == "g1"   # gene column appended

    # minimap2 index built; hisat2 skipped
    assert res["ref_mm2_index"] and os.path.isfile(res["ref_mm2_index"])
    assert res["hisat2_index"] is None


def test_cli_init_ref(tmp_path):
    fa, gtf = _make_genome(tmp_path)
    rc = main(["init-ref", "-f", fa, "--gtf", gtf, "-o", str(tmp_path)])
    assert rc == 0
    assert os.path.isfile(str(tmp_path / "genome.len.txt"))
    assert os.path.isfile(str(tmp_path / "genome.mm2"))


def _ensembl_gff3(gtf, out):
    """The GTF as an Ensembl-style GFF3: gene/mRNA/exon lines whose IDs carry the gene:/transcript: type prefixes."""
    tx = {}
    for line in open(gtf):
        f = line.rstrip("\n").split("\t")
        gene, t = re.search(r'gene_id "([^"]+)"', f[8]).group(1), re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        tx.setdefault((gene, t), {"chrom": f[0], "strand": f[6], "exons": []})["exons"].append((int(f[3]), int(f[4])))
    genes = {}
    for (gene, t), v in tx.items():
        g = genes.setdefault(gene, {"chrom": v["chrom"], "strand": v["strand"], "start": 10**12, "end": 0, "tx": []})
        g["start"], g["end"] = min(g["start"], v["exons"][0][0]), max(g["end"], v["exons"][-1][1])
        g["tx"].append((t, v))
    with open(out, "w") as fh:
        fh.write("##gff-version 3\n")
        for gene, g in genes.items():
            fh.write(f"{g['chrom']}\tens\tgene\t{g['start']}\t{g['end']}\t.\t{g['strand']}\t.\tID=gene:{gene};biotype=protein_coding\n")
            for t, v in g["tx"]:
                fh.write(f"{v['chrom']}\tens\tmRNA\t{v['exons'][0][0]}\t{v['exons'][-1][1]}\t.\t{v['strand']}\t.\t"
                         f"ID=transcript:{t};Parent=gene:{gene}\n")
                for a, b in v["exons"]:
                    fh.write(f"{v['chrom']}\tens\texon\t{a}\t{b}\t.\t{v['strand']}\t.\tParent=transcript:{t}\n")


@pytest.mark.skipif(shutil.which("gffread") is None, reason="needs gffread (core iflas env)")
def test_init_ref_from_ensembl_gff3_matches_the_gtf_route(tmp_path):
    """`iflas init-ref --gff` converts with gffread and drops the Ensembl gene:/transcript: prefixes: the derived
    gene models equal those of the GTF the GFF3 was written from (gffread was missing from the env: the route failed)."""
    gff = tmp_path / "genes.gff3"
    _ensembl_gff3(os.path.join(V2, "genes.gtf"), gff)
    shutil.copy(os.path.join(V2, "genome.fa"), tmp_path / "genome.fa")
    from_gtf = init_ref(str(tmp_path / "genome.fa"), gtf=os.path.join(V2, "genes.gtf"), outdir=str(tmp_path / "gtf"), build_minimap2=False)
    rc = main(["init-ref", "-f", str(tmp_path / "genome.fa"), "--gff", str(gff), "-o", str(tmp_path / "gff")])
    assert rc == 0
    def models(gpe):
        return sorted(tuple(f[i] for i in (0, 1, 2, 3, 4, 7, 8, 9, 11)) for f in (l.split("\t") for l in open(gpe).read().splitlines()))
    assert models(tmp_path / "gff" / "genome.gpe") == models(from_gtf["ref_gpe"])
    assert ids_with_colons(str(tmp_path / "gff" / "genome.gtf")) == []


@pytest.mark.skipif(shutil.which("gffread") is None or shutil.which("hisat2-build") is None, reason="needs the core iflas env")
def test_init_ref_decompresses_gzipped_inputs(tmp_path):
    """BUG-194: download sites ship .fa.gz/.gff3.gz. A bgzip genome with a gzipped GFF3, and a plain-gzip genome with
    a gzipped GTF, give the same gene models as the plain files, from uncompressed copies in the output dir (before:
    gffread read the .gff3.gz as empty, the .gtf.gz crashed with UnicodeDecodeError, samtools faidx refused plain
    gzip and hisat2-build any compressed FASTA)."""
    import gzip
    import subprocess
    gff = tmp_path / "genes.gff3"
    _ensembl_gff3(os.path.join(V2, "genes.gtf"), gff)
    genome = os.path.join(V2, "genome.fa")
    subprocess.run(f"bgzip -c {genome} > {tmp_path}/genome.fa.gz && gzip -c {gff} > {tmp_path}/genes.gff3.gz", shell=True, check=True)
    with open(genome, "rb") as src, gzip.open(tmp_path / "plain_gzip.fa.gz", "wb") as dst:
        dst.write(src.read())
    with open(os.path.join(V2, "genes.gtf"), "rb") as src, gzip.open(tmp_path / "genes.gtf.gz", "wb") as dst:
        dst.write(src.read())
    plain = init_ref(genome, gtf=os.path.join(V2, "genes.gtf"), outdir=str(tmp_path / "plain"), build_minimap2=False)
    from iflas.runtime import ToolRunner

    class Recording(ToolRunner):
        commands = []

        def run(self, cmd, **kwargs):
            self.commands.append(cmd)
            return super().run(cmd, **kwargs)
    bgz = init_ref(str(tmp_path / "genome.fa.gz"), gff=str(tmp_path / "genes.gff3.gz"), outdir=str(tmp_path / "bgz"),
                   build_minimap2=False, build_hisat2=True, threads=2, runner=Recording())
    assert ["hisat2-build", "-q", "-p", "2"] == next(c for c in Recording.commands if c[0] == "hisat2-build")[:4]     # init-ref -t
    gz = init_ref(str(tmp_path / "plain_gzip.fa.gz"), gtf=str(tmp_path / "genes.gtf.gz"), outdir=str(tmp_path / "gz"),
                  build_minimap2=False)

    def models(gpe):
        return sorted(tuple(f[i] for i in (0, 1, 2, 3, 4, 7, 8, 9, 11)) for f in (l.split("\t") for l in open(gpe).read().splitlines()))
    assert models(bgz["ref_gpe"]) == models(gz["ref_gpe"]) == models(plain["ref_gpe"])
    assert bgz["ref_genome"] == str(tmp_path / "bgz" / "genome.fa") and gz["ref_gtf"] == str(tmp_path / "gz" / "genes.gtf")
    assert open(bgz["ref_genome"]).read() == open(gz["ref_genome"]).read() == open(genome).read()
    assert open(bgz["ref_size"]).read() == open(plain["ref_size"]).read() and os.path.isfile(bgz["hisat2_index"] + ".1.ht2")


def test_ids_with_colons_lists_ids_iflas_cannot_split(tmp_path):
    gtf = tmp_path / "g.gtf"
    gtf.write_text('c\ts\texon\t1\t9\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
                   'c\ts\texon\t1\t9\t.\t+\t.\tgene_id "LOC:2"; transcript_id "rna:2";\n')
    assert ids_with_colons(str(gtf)) == ["LOC:2", "rna:2"]
