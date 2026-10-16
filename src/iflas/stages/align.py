"""Alignment stage core — TGS via minimap2, NGS via hisat2 (ported from legacy mapping.py).

Faithful alignment parameters from the published pipeline. Includes the NGS reassembly +
junction extraction (legacy `srMapping` + `getJuncFromRegtools`) that produces
`reassembly/tmp.bam` + `junctions.bed` (consumed by pu_filter feature extraction and the
find_as NGS confident-filtering). The iFLAS-specific junction-combination support filtering
(legacy `mappingFilterAndAddTags`) is available below as `filter_by_junction` /
`mapping_filter_and_add_tags` and is covered by command-contract tests.
"""

from __future__ import annotations

import glob
import os
import tempfile

from ..io.files import removeFiles
from ..io.records import BedFile
from ..runtime import ToolRunner, quote as q


BAI_MAX_LENGTH = 2 ** 29 - 1             # 536,870,911 bp: the longest reference a .bai index can address
MINIMAP2_SPLIT_BASES = 4_000_000_000     # minimap2 cuts its index every -I bases (8G in 2.31, 4G in older releases)


def _ensure_parent(path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)


def longest_reference(bam) -> int:
    """Length of the longest @SQ of a BAM header (0 when it cannot be read yet)."""
    import pysam

    try:
        with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as f:
            return max(f.lengths, default=0)
    except (OSError, ValueError):
        return 0


def index_bam(bam, *, runner, cwd=None):
    """``samtools index``: a .bai, or a .csi when a reference is longer than a .bai can address (wheat, barley, ...;
    BUG-189)."""
    args = ["samtools", "index"] + (["-c"] if longest_reference(bam) > BAI_MAX_LENGTH else []) + [bam]
    return runner.run(args, **({"cwd": cwd} if cwd else {}))


def minimap2_index_batch(genome_bases: int) -> str | None:
    """``-I`` that keeps minimap2's index in one part for a genome above MINIMAP2_SPLIT_BASES: a multi-part index writes
    SAM without @SQ lines and reports each read once per part, and --split-prefix, minimap2's remedy, refuses --MD
    (BUG-190). None for smaller genomes (command unchanged)."""
    if genome_bases <= MINIMAP2_SPLIT_BASES:
        return None
    return f"{-(-genome_bases // 1_000_000_000) + 1}G"


def genome_bases(sizes_file) -> int:
    """Total length from a chrom<TAB>length file (ref_size or a .fai); 0 when there is none."""
    if not sizes_file or not os.path.isfile(sizes_file):
        return 0
    with open(sizes_file) as f:
        return sum(int(line.split("\t")[1]) for line in f if line.strip())


def check_single_part_index(sam_file) -> None:
    """Stop on minimap2 output from a multi-part index (no @SQ header lines) instead of filtering corrupt alignments."""
    with open(sam_file) as f:
        for line in f:
            if line.startswith("@SQ"):
                return
            if not line.startswith("@"):
                raise ValueError(f"{sam_file}: no @SQ header lines: the minimap2 index is split into parts (genome larger "
                                 "than its -I). Rebuild it with `iflas init-ref` or `minimap2 -I <genome size, e.g. 20G> -d`")


def minimap2_tgs(target, flnc, out_sam, *, platform="pacbio", max_intron=50000,
                 threads=4, index_batch=None, runner=None):
    """Splice-align full-length TGS reads to `target` (genome FASTA or .mm2 index) -> SAM.

    pacbio:   minimap2 -ax splice:hq -G <max_intron> -uf --secondary=no --MD
    nanopore: minimap2 -ax splice    -G <max_intron> -k14 -uf --secondary=no --MD
    The SAM feeds ``mapping_filter_and_add_tags`` (legacy ``lrMapping``), which writes the sorted BAM.
    ``index_batch`` (``minimap2_index_batch``) keeps the index of a large genome in one part.
    """
    runner = runner or ToolRunner()
    _ensure_parent(out_sam)
    if platform == "pacbio":
        preset, extra = "splice:hq", ""
    else:
        preset, extra = "splice", " -k14"
    if index_batch:
        extra += f" -I {index_batch}"
    runner.run(f"minimap2 -ax {preset} -G {max_intron} -uf{extra} --secondary=no --MD "
               f"{q(target)} {q(flnc)} -t {threads} > {q(out_sam)}")
    return out_sam


def flnc_to_fasta(flnc, out_fa, *, runner=None):
    """``seqkit fq2fa -w 0`` -> legacy ``mapping/flnc.processed.fa`` (collapse input; FASTA passes through)."""
    runner = runner or ToolRunner()
    _ensure_parent(out_fa)
    runner.run(f"seqkit fq2fa -w 0 {q(flnc)} > {q(out_fa)}")
    return out_fa


def _without_spaces(paths, tmp):
    """``paths`` (comma-separated files) with every file whose path holds whitespace replaced by a symlink in ``tmp``:
    the hisat2 wrapper re-splits its arguments on spaces, so quoting alone does not help."""
    out = []
    for path in paths.split(","):
        if any(c.isspace() for c in path):
            link = os.path.join(tmp, f"{len(os.listdir(tmp))}_{os.path.basename(path).replace(' ', '_')}")
            os.symlink(os.path.abspath(path), link)
            path = link
        out.append(path)
    return ",".join(out)


def hisat2_ngs(index, out_bam, *, left=None, right=None, single=None, threads=4,
               max_intron=50000, runner=None):
    """Align NGS reads with HISAT2 (--dta) -> sorted, indexed BAM. ``--reorder`` keeps the input order of the reads, so
    the BAM -- and the StringTie assembly and regtools junction names built from it -- no longer change with thread
    scheduling (BUG-163; legacy wrote them in completion order)."""
    runner = runner or ToolRunner()
    _ensure_parent(out_bam)
    if not ((left and right) or single or left):
        raise ValueError("hisat2_ngs needs paired (left+right) or single reads")
    # the links must sit in a path without whitespace: the configured tmp_dir (BUG-184) unless it has some
    link_dir = tempfile.gettempdir() if not any(c.isspace() for c in tempfile.gettempdir()) else "/tmp"
    with tempfile.TemporaryDirectory(prefix="iflas_hisat2_", dir=link_dir) as tmp:
        if any(c.isspace() for c in str(index)):
            for f in glob.glob(glob.escape(str(index)) + ".*.ht2*"):
                os.symlink(os.path.abspath(f), os.path.join(tmp, "index" + f[len(str(index)):]))
            index = os.path.join(tmp, "index")
        if left and right:
            reads = f"-1 {q(_without_spaces(left, tmp))} -2 {q(_without_spaces(right, tmp))}"
        else:
            reads = f"-U {q(_without_spaces(single or left, tmp))}"
        runner.run(f"hisat2 -x {q(index)} {reads} --dta -p {threads} --reorder --max-intronlen {max_intron} "
                   f"| samtools sort -@ {threads} -o {q(out_bam)} -")
    index_bam(out_bam, runner=runner)
    return out_bam


def ngs_reassembly_junctions(hisat2_index, out_dir, *, repeats, threads=4,
                             max_intron=50000, filter_count=10, runner=None):
    """NGS reassembly + junction extraction (ported from srMapping + getJuncFromRegtools).

    For each repeat (dict ``{"left","right"}`` or ``{"single"}``): HISAT2 -> sorted BAM,
    then per-repeat ``stringtie``. Repeats are then merged (``stringtie --merge``) and the
    BAMs concatenated+sorted into ``reassembly/tmp.bam``; ``regtools junctions extract``
    yields junctions, kept when read support (col5) > ``filter_count`` (legacy
    ``awk '$5>count'``). Returns the tmp_bam / junctions_bed / merged_gtf paths.
    """
    runner = runner or ToolRunner()
    align_dir = os.path.join(out_dir, "alignment")
    reasm_dir = os.path.join(out_dir, "reassembly")
    os.makedirs(align_dir, exist_ok=True)
    os.makedirs(reasm_dir, exist_ok=True)

    bam_list, gtf_list = [], []
    for i, rep in enumerate(repeats):
        name = "repeat" + str(i)
        bam = os.path.join(align_dir, name + ".sorted.bam")
        hisat2_ngs(hisat2_index, bam, left=rep.get("left"), right=rep.get("right"),
                   single=rep.get("single"), threads=threads, max_intron=max_intron, runner=runner)
        gtf = os.path.join(align_dir, name + ".gtf")
        runner.run(["stringtie", bam, "-o", gtf, "-p", str(threads)])
        bam_list.append(bam)
        gtf_list.append(gtf)

    merged_gtf = os.path.join(reasm_dir, "stringtie_merged.gtf")
    runner.run(["stringtie", "--merge", "-p", str(threads), "-o", merged_gtf, *gtf_list])

    tmp_bam = os.path.join(reasm_dir, "tmp.bam")
    runner.run(f"samtools cat {' '.join(map(q, bam_list))} | samtools sort -@ {threads} > {q(tmp_bam)}")
    index_bam(tmp_bam, runner=runner)

    raw_bed = os.path.join(reasm_dir, "junctions.raw.bed")
    junctions_bed = os.path.join(reasm_dir, "junctions.bed")
    # regtools 1.0 CLI delta: strandness is now a string {XS,RF,FR} (was numeric 0/1/2);
    # legacy `-s 0` (unstranded, use aligner tags) -> `-s XS` (hisat2 --dta emits XS).
    runner.run(["regtools", "junctions", "extract", "-a", "5", "-m", "50", "-M", str(max_intron),
                tmp_bam, "-s", "XS", "-o", raw_bed])
    with open(raw_bed) as fin, open(junctions_bed, "w") as fout:
        for line in fin:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 5 and int(cols[4]) > filter_count:
                fout.write(line)
    return {"tmp_bam": tmp_bam, "junctions_bed": junctions_bed, "merged_gtf": merged_gtf}


def _supported(read, junc_count, comb_reads, ref_juncs, ngs_juncs, max_intron_len, junc_comb_support_n):
    """A multi-exon read passes with every intron within the maximum, its junction chain in at least
    ``junc_comb_support_n`` reads, and either >= 4 introns or >= 1/3 consensus junctions (reference, NGS or >= 3 FLNC).
    BUG-193 (adopted 2026-10-16): the published code compared the total intron length, which dropped every full-length
    read of a long multi-intron gene (162 maize B73 transcripts)."""
    if max(end - start for start, end in read.introns) > max_intron_len or len(comb_reads[read.juncChain]) < junc_comb_support_n:
        return False
    if len(read.introns) >= 4:
        return True
    consensus = sum(1 for junc in read.introns if junc in ref_juncs or junc in ngs_juncs or junc_count[junc] >= 3)
    return consensus / float(len(read.introns)) >= 1.0 / 3


def filter_by_junction(bed_file=None, out_prefix="flnc", ref_bed_file=None,
                       junc_bed_file=None, max_intron_len=50000,
                       junc_comb_support_n=2, work_dir=None):
    """Filter FLNC alignments by legacy junction-combination support.

    Port of legacy ``mapping.py::filterByJunc`` (intron length per intron since BUG-193).  Multi-exon reads are retained
    when their junction chain is sufficiently supported and either enough junctions are
    consensus (reference, NGS, or >=3 FLNC support) or the read has >=4 introns.
    Single-exon reads pass through unchanged.
    """
    cwd = os.getcwd()
    if work_dir is not None:
        os.makedirs(work_dir, exist_ok=True)
        os.chdir(work_dir)
    try:
        ref_bed = BedFile(ref_bed_file, type="bed12+").reads
        flnc_bed = BedFile(bed_file, type="bed12+").reads
        ngs_juncs = set()
        if junc_bed_file:
            ngs_juncs = {junc for r in BedFile(junc_bed_file, type="bed12").reads.values() for junc in r.introns}
        ref_juncs = {junc for r in ref_bed.values() for junc in r.introns}
        junc_count, comb_reads = {}, {}
        for name, read in flnc_bed.items():
            if read.introns:
                for junc in read.introns:
                    junc_count[junc] = junc_count.get(junc, 0) + 1
                comb_reads.setdefault(read.juncChain, []).append(name)

        out_file = "{}.filterReads.lst".format(out_prefix)
        with open(out_file, "w") as filter_reads:
            for name, read in flnc_bed.items():
                if not read.introns or _supported(read, junc_count, comb_reads, ref_juncs, ngs_juncs, max_intron_len,
                                                  junc_comb_support_n):
                    print(name, file=filter_reads)
        return os.path.join(os.getcwd(), out_file)
    finally:
        if work_dir is not None:
            os.chdir(cwd)


def mapping_filter_and_add_tags(sam_file=None, out_prefix="flnc", max_length=50000,
                                ref_bed_file=None, junc_bed_file=None, threads=10,
                                junc_comb_sup=2, work_dir=None, runner=None):
    """Port legacy ``mappingFilterAndAddTags`` with ToolRunner-managed commands."""
    runner = runner or ToolRunner()
    work_dir = os.path.abspath(work_dir or os.getcwd())
    os.makedirs(work_dir, exist_ok=True)

    script_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
    util_dir = os.path.join(script_dir, "utils")
    sam_file = os.path.abspath(sam_file)
    if os.path.isfile(sam_file):
        check_single_part_index(sam_file)
    ref_bed_file = os.path.abspath(ref_bed_file) if ref_bed_file else None
    junc_bed_file = os.path.abspath(junc_bed_file) if junc_bed_file else None

    runner.run(["bash", "-c", f"bamToBed -i <(samtools view -bS {q(sam_file)}) -bed12 > tmp.bed12"],
               cwd=work_dir)
    filtered_reads = filter_by_junction(
        bed_file=os.path.join(work_dir, "tmp.bed12"),
        out_prefix=out_prefix,
        ref_bed_file=ref_bed_file,
        junc_bed_file=junc_bed_file,
        max_intron_len=max_length,
        junc_comb_support_n=junc_comb_sup,
        work_dir=work_dir,
    )
    runner.run(
        "(samtools view -H {s}; {u}/filter.pl -o {r} {s} -m i) > tmp.sam".format(s=q(sam_file), u=q(util_dir), r=q(filtered_reads)),
        cwd=work_dir, shell=True,
    )
    runner.run(
        "(samtools view -H tmp.sam; samtools view -f 16 -F 4079 tmp.sam; "
        "samtools view -f 0 -F 4095 tmp.sam) | "
        "samtools sort -@ {} --output-fmt SAM > {}.mm2.sam".format(threads, out_prefix),
        cwd=work_dir, shell=True,
    )
    runner.run(
        "{}/samAddTag.pl --checkHardClip --coverage --identity --u unmmaped.sam {}.mm2.sam "
        "2> lengthInconsistent.sam | {}/sam2bed.pl -t CV,ID > {}.addCVandID.bed12+ "
        "2>/dev/null".format(q(util_dir), out_prefix, q(util_dir), out_prefix),
        cwd=work_dir, shell=True,
    )
    runner.run(
        "samtools view -bS -@ {} {}.mm2.sam > {}.mm2.sorted.bam".format(
            threads, out_prefix, out_prefix,
        ),
        cwd=work_dir, shell=True,
    )
    index_bam(os.path.join(work_dir, "{}.mm2.sorted.bam".format(out_prefix)), runner=runner, cwd=work_dir)
    removeFiles(work_dir, ["tmp.bed12", "tmp.sam"])
    return {
        "filtered_reads": filtered_reads,
        "sam": os.path.join(work_dir, "{}.mm2.sam".format(out_prefix)),
        "bed12plus": os.path.join(work_dir, "{}.addCVandID.bed12+".format(out_prefix)),
        "bam": os.path.join(work_dir, "{}.mm2.sorted.bam".format(out_prefix)),
    }
