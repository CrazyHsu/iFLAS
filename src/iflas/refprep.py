"""Reference preparation: derive aux files from genome FASTA + GTF.

Replaces the long manual list of absolute-path reference files (ref_size / ref_gpe /
ref_bed / ref_mm2_index / hisat2_index) — the user supplies only genome + GTF.
Also fixes the legacy initRefSetting size derivation (BUG-110: opened the size file
read-only and wrote the sequence instead of its length).
"""

from __future__ import annotations

import gzip
import os
import re
import shutil

from .io.genepred import GenePredObj
from .runtime import ToolRunner, get_logger


def _stem(path: str) -> str:
    return os.path.splitext(os.path.basename(path))[0]


def is_gzipped(path: str) -> bool:
    """True for gzip and bgzip files (both start with the gzip magic bytes)."""
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def plain_copy(path: str, outdir: str, logger) -> str:
    """``path``, or for a gzip/bgzip file its decompressed copy in ``outdir``: gffread reads a compressed GFF3 as empty,
    samtools faidx refuses plain gzip and hisat2-build any compressed FASTA (BUG-194). A copy at least as new as
    ``path`` is reused."""
    if not os.path.isfile(path) or not is_gzipped(path):
        return path
    name = os.path.basename(path)
    target = os.path.join(outdir, name[:-3] if name.endswith(".gz") else name + ".plain")
    if not (os.path.isfile(target) and os.path.getmtime(target) >= os.path.getmtime(path)):
        logger.info("decompressing %s -> %s", path, target)
        with gzip.open(path, "rb") as src, open(target, "wb") as dst:
            shutil.copyfileobj(src, dst, 1 << 20)
    return target


def write_chrom_sizes(fai_path: str, out_path: str) -> None:
    """Write a 2-column chrom<TAB>length file from a samtools .fai index."""
    with open(fai_path) as fai, open(out_path, "w") as out:
        for line in fai:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                out.write(f"{parts[0]}\t{parts[1]}\n")


_ENSEMBL_ID_PREFIX = re.compile(r'\b(gene_id|transcript_id) "(?:gene|transcript):')
_GTF_ID = re.compile(r'\b(?:gene_id|transcript_id) "([^"]*)"')


def strip_ensembl_id_prefixes(gtf_path: str) -> int:
    """Drop the ``gene:``/``transcript:`` type prefixes gffread keeps from Ensembl GFF3 IDs (Ensembl's own GTFs carry the
    bare IDs, and iFLAS uses ':' as a separator in event names). Rewrites ``gtf_path``; returns the lines changed."""
    with open(gtf_path) as f:
        lines = f.readlines()
    stripped = [_ENSEMBL_ID_PREFIX.sub(r'\1 "', line) for line in lines]
    with open(gtf_path, "w") as f:
        f.writelines(stripped)
    return sum(a != b for a, b in zip(lines, stripped))


def ids_with_colons(gtf_path: str, limit: int = 5) -> list[str]:
    """Up to ``limit`` gene/transcript IDs containing ':', which iFLAS splits event and junction names on."""
    found: list[str] = []
    with open(gtf_path) as f:
        for line in f:
            for value in _GTF_ID.findall(line):
                if ":" in value and value not in found:
                    found.append(value)
                    if len(found) == limit:
                        return found
    return found


def gpe_to_bed(gpe_path: str, out_path: str, gene: bool = True) -> None:
    """Convert a gpe file to BED12(+gene) via the ported GenePredObj (no perl gpe2bed)."""
    obj = GenePredObj(gpe_path, bincolumn=False)
    with open(out_path, "w") as out:
        for chrom in obj.genePredDict:
            for strand in obj.genePredDict[chrom]:
                for gp in obj.genePredDict[chrom][strand]:
                    out.write(str(gp[2].go_to_bed(gene=gene)) + "\n")


def init_ref(genome: str, gtf: str | None = None, gff: str | None = None,
             outdir: str | None = None, *, build_minimap2: bool = True,
             build_hisat2: bool = False, threads: int = 1, runner: ToolRunner | None = None, logger=None) -> dict:
    """Derive reference aux files; returns a dict of the reference fields for config.yaml."""
    logger = logger or get_logger()
    runner = runner or ToolRunner(logger=logger)

    genome = os.path.abspath(genome)
    if not os.path.isfile(genome):
        raise FileNotFoundError(genome)
    outdir = os.path.abspath(outdir) if outdir else os.path.dirname(genome)
    os.makedirs(outdir, exist_ok=True)
    genome = plain_copy(genome, outdir, logger)
    base = os.path.join(outdir, _stem(genome))

    if gtf is None and gff is not None:
        gtf = base + ".gtf"
        runner.run(["gffread", "-T", "-o", gtf, plain_copy(os.path.abspath(gff), outdir, logger)])
        if os.path.isfile(gtf):
            strip_ensembl_id_prefixes(gtf)
    if gtf is None:
        raise ValueError("init_ref needs a GTF (or a GFF to convert via gffread)")
    gtf = plain_copy(os.path.abspath(gtf), outdir, logger)
    odd = ids_with_colons(gtf) if os.path.isfile(gtf) else []
    if odd:
        logger.warning("gene/transcript IDs contain ':' (e.g. %s); iFLAS joins event and junction names with ':', "
                       "so rename these IDs in the annotation", ", ".join(odd))

    # chromosome sizes via samtools faidx
    runner.run(["samtools", "faidx", genome])
    size = base + ".len.txt"
    write_chrom_sizes(genome + ".fai", size)

    # gene predictions (gpe) via UCSC gtfToGenePred
    gpe = base + ".gpe"
    runner.run(["gtfToGenePred", gtf, gpe, "-genePredExt"])

    # BED12(+gene) via the ported GenePredObj (replaces utils/gpe2bed.pl; verify vs baseline)
    bed = base + ".bed"
    gpe_to_bed(gpe, bed, gene=True)

    result = {"ref_genome": genome, "ref_gtf": gtf, "ref_size": size,
              "ref_gpe": gpe, "ref_bed": bed, "ref_mm2_index": None, "hisat2_index": None}

    if build_minimap2:
        from .stages.align import genome_bases, minimap2_index_batch

        mm2 = base + ".mm2"
        batch = minimap2_index_batch(genome_bases(size))
        runner.run(["minimap2"] + (["-I", batch] if batch else []) + ["-d", mm2, genome])
        result["ref_mm2_index"] = mm2
    if build_hisat2:
        prefix = base + "_hisat2"
        runner.run(["hisat2-build", "-q"] + (["-p", str(threads)] if threads > 1 else []) + [genome, prefix])
        result["hisat2_index"] = prefix

    logger.info("init-ref done: %s", base)
    return result
