"""Collapse stage — redundant-isoform collapse via cDNA_Cupcake (ported from legacy collapse.py).

cupcake pins Python 3.10, so it lives in the isolated `iflas-cupcake` env and is invoked via
`conda run -n iflas-cupcake collapse_isoforms_by_sam.py`. cupcake 29 accepts a sorted BAM (-b)
directly, so no SAM conversion is needed. Parameters are faithful to the published pipeline.
"""

from __future__ import annotations

import os

from ..runtime import ToolError, ToolRunner


def collapse_isoforms(flnc_fa, sorted_bam, out_prefix, *, max_5_diff=1000, max_3_diff=500,
                      fl_coverage=2, min_coverage=0.9, min_identity=0.9, max_fuzzy_junction=5,
                      dun_merge_5_shorter=True, conda_env="iflas-cupcake", runner=None):
    """Collapse aligned full-length reads into non-redundant isoforms.

    Returns the path to the collapsed GFF (``<out_prefix>.collapsed.gff``).
    """
    runner = runner or ToolRunner()
    cmd = [
        "collapse_isoforms_by_sam.py",
        "--input", flnc_fa,
        "-b", sorted_bam,
        "-o", out_prefix,
        "-c", str(min_coverage),
        "-i", str(min_identity),
        "--max_fuzzy_junction", str(max_fuzzy_junction),
        "--max_5_diff", str(max_5_diff),
        "--max_3_diff", str(max_3_diff),
        "--flnc_coverage", str(fl_coverage),
    ]
    if dun_merge_5_shorter:
        cmd.append("--dun-merge-5-shorter")
    # cupcake 29 exits non-zero on the empty .collapsed.bad.gff GFF3 sanity check even when the
    # real (.good.gff) output is fine, so don't trust the return code — verify the output file.
    proc = runner.run(cmd, conda_env=conda_env, check=False)
    for cand in (out_prefix + ".collapsed.good.gff", out_prefix + ".collapsed.gff"):
        if os.path.isfile(cand) and os.path.getsize(cand) > 0:
            return cand
    raise ToolError("collapse_isoforms_by_sam.py", proc.returncode,
                    (getattr(proc, "stderr", "") or "no collapsed GFF produced"))
