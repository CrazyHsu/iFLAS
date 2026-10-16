"""iFLAS command-line interface."""

from __future__ import annotations

import argparse
from pathlib import Path
import os
import sys

from . import __version__


def _count_nonempty_lines(path, *, missing_ok=False) -> int:
    if missing_ok and not os.path.isfile(path):
        return 0
    with open(path) as handle:
        return sum(1 for line in handle if line.strip())


_REPORT_PLOTS = (      # (flag, done message, plotter in report.plots, message when nothing was found, takes --sample)
    ("basic", "Generated basic report plots", "plot_basic_report", "no FLNC/mapping inputs found", True),
    ("asp", "Generated AS pattern report plots", "plot_as_pattern_report", "no AS characterization inputs found", True),
    ("gene_struc", "Generated target gene structure report plots", "plot_target_gene_structure_report",
     "no target gene structure inputs found", True),
    ("novel_hq_as", "Generated report plots", "plot_novel_hq_as_report", "no isoformScore.txt files found", True),
    ("diff", "Generated differential AS report plots", "plot_diff_as_report", "no das/*.sigDiffAS directories found", False),
    ("asas", "Generated allele-specific AS report plots", "plot_asas_report", "no ASAS report inputs found", True),
    ("palen", "Generated poly(A) tail report plots", "plot_palen_report", "no palenAS report inputs found", True),
)


def _cmd_report(args) -> int:
    wanted = [row for row in _REPORT_PLOTS if args.all or getattr(args, row[0])]
    run_html = args.html or args.all
    if not wanted and not run_html:
        sys.stderr.write("iflas report: no output selected; pass -all or one of -basic/-asp/-geneStruc/--novel-hq-as/"
                         "-diff/-asas/-palen/-html.\n")
        return 2
    if not args.config or not args.samples:
        sys.stderr.write("iflas report: provide -c/--config and -s/--samples for report outputs.\n")
        return 2

    from .report import plots as report_plots

    try:
        cfg = _load_config(args)
        for flag, done, plotter, empty, per_sample in wanted:
            fn = getattr(report_plots, plotter)
            plots = fn(cfg, sample_names=args.sample) if per_sample else fn(cfg)
            print(done)
            for plot in plots["plots"]:
                print(f"  plot: {plot}")
            if not plots["plots"]:
                print(f"  {empty}")
        if run_html:
            from .report.html import generate_html_report

            out = generate_html_report(cfg, sample_names=args.sample)
            print("Wrote HTML report")
            print(f"  html: {out['main_page']}")
            print(f"  json: {out['json']}")
    except (FileNotFoundError, ValueError) as exc:
        sys.stderr.write(f"iflas report: {exc}\n")
        return 2
    return 0


def _load_config(args):
    """The run's config; a set ``dirs.tmp_dir`` receives the temporary files from here on (BUG-184)."""
    from .config import load_config
    from .runtime import use_tmp_dir

    cfg = load_config(args.config, args.samples)
    if cfg.dirs.tmp_dir:
        use_tmp_dir(cfg.dirs.tmp_dir)
    return cfg


def _cmd_migrate_config(args) -> int:
    from .config import migrate_config_file
    config, samples, skipped = migrate_config_file(args.input, args.out, args.samples)
    print(f"Wrote {args.out} (reference entries: {len(config.get('reference', []))})")
    print(f"Wrote {args.samples} ({len(samples)} samples)")
    if skipped:
        print(f"Skipped sections: {', '.join(skipped)}")
    return 0


def _cmd_init_ref(args) -> int:
    from .refprep import init_ref
    res = init_ref(args.fasta, gtf=args.gtf, gff=args.gff, outdir=args.outdir,
                   build_minimap2=not args.no_minimap2, build_hisat2=args.hisat2, threads=args.threads)
    print("Derived reference files (paste under 'reference:' in config.yaml):")
    for key, val in res.items():
        if val:
            print(f"  {key}: {val}")
    return 0


def _cmd_find_as(args) -> int:
    from .stages.find_as import find_as_lr, find_as_sample

    outdir = args.outdir or "."
    if args.base_dir:
        out = find_as_sample(
            args.base_dir, args.ref_gpe, out_dir=args.outdir,
            isoform_bed12plus=args.isoform_bed, group_file=args.group,
            ngs_junctions=args.ngs_junctions)
    elif args.isoform_bed and args.group:
        out = find_as_lr(args.isoform_bed, args.group, args.ref_gpe, outdir,
                         ngs_junctions=args.ngs_junctions)
    else:
        sys.stderr.write("iflas find_as: provide --base-dir or both --isoform-bed and --group.\n")
        return 2

    print("Wrote AS events (LR):")
    for as_type in ("IR", "SE", "A5SS", "A3SS"):
        path = out[as_type]
        n = _count_nonempty_lines(path)
        print(f"  {as_type}: {path} ({n} events)")
    print("Wrote confident AS events:")
    for as_type in ("IR", "SE", "A5SS", "A3SS"):
        path = out["confident"][as_type]
        n = _count_nonempty_lines(path)
        print(f"  {as_type}: {path} ({n} events)")
    return 0


def _cmd_find_pa(args) -> int:
    from .stages.find_pa import find_pa_sample, getPaCluster

    if args.base_dir:
        if not (args.ref_genome and args.ref_size):
            sys.stderr.write("iflas find_pa: --base-dir mode needs --ref-genome and --ref-size.\n")
            return 2
        out = find_pa_sample(args.base_dir, args.ref_genome, args.ref_size, filter_pa_by_rpkm=args.pa_rpkm,
                             filter_pa_by_count=args.pa_support, confident_pa=args.confident_pa,
                             threads=args.threads or 1)
        print(f"Wrote {out['pa_bed6plus']} ({_count_nonempty_lines(out['pa_bed6plus'])} PA sites); motif plots in {out['motif_dir']}")
        return 0
    if not (args.reads_bed and args.group and args.out):
        sys.stderr.write("iflas find_pa: provide --base-dir (whole analysis) or --reads-bed, --group and -o (clustering only).\n")
        return 2
    getPaCluster(
        readsBed=args.reads_bed,
        tofuGroup=args.group,
        filterByCount=args.pa_support,
        paClusterOut=args.out,
        paDist=args.pa_dist,
        windowSize=args.window_size,
        rpkmPAC=args.pa_rpkm,
        confidentPa=args.confident_pa,
    )
    n = _count_nonempty_lines(args.out, missing_ok=True)
    print(f"Wrote {args.out} ({n} PA clusters)")
    return 0


def _split_cli_paths(values):
    paths = []
    for value in values or []:
        paths.extend(part.strip() for part in str(value).split(",") if part.strip())
    return paths


def _cmd_palen_as(args) -> int:
    from .stages.palen_as import palen_as_sample

    try:
        out = palen_as_sample(
            args.base_dir,
            _split_cli_paths(args.palen),
            fast5_dir=args.fast5_dir,
            ref_genome=args.ref_genome,
            sequencing_summary=args.sequencing_summary,
            threads=args.threads or 1,
            out_dir=args.outdir,
            isoform_file=args.isoform_bed,
            group_file=args.group,
            reads_file=args.reads_bed,
            as_dir=args.as_dir,
            pa_cluster=args.pa_cluster,
            filter_by_count=args.pa_support,
            merge_by_junc=not args.no_merge_by_junc,
            confident_pac_files=_split_cli_paths(args.confident_pac),
            ref_bed=args.ref_bed,
        )
    except (FileNotFoundError, ValueError, KeyError) as exc:
        sys.stderr.write(f"iflas palen_as: {exc}\n")
        return 2

    print("Wrote palen_as outputs:")
    print(f"  AS dir: {out['as_dir']}")
    for as_type, path in sorted(out["as_sig"].items()):
        n = _count_nonempty_lines(path, missing_ok=True)
        print(f"  {as_type}: {path} ({n} records)")
    if out.get("apa"):
        n = _count_nonempty_lines(out["apa"], missing_ok=True)
        print(f"  APA: {out['apa']} ({n} records)")
    return 0


def _cmd_asas(args) -> int:
    from .stages.asas import asas_sample

    try:
        out = asas_sample(
            args.base_dir,
            out_dir=args.outdir,
            isoform_file=args.isoform_bed,
            group_file=args.group,
            as_dir=args.as_dir,
            loci_dir=args.loci_dir,
            filter_by_count=args.filter_by_count,
            merge_by_junc=args.merge_by_junc,
            ref_genome=args.ref_genome,
            flnc=args.flnc,
            ngs_confirm=not args.no_ngs_confirm,
            collapsed_gff=args.collapsed_gff,
            min_coverage=args.min_coverage,
            threads=args.threads or 1,
        )
    except (FileNotFoundError, KeyError, ValueError) as exc:
        sys.stderr.write(f"iflas asas: {exc}\n")
        return 2

    n = _count_nonempty_lines(out["partial"], missing_ok=True)
    print("Wrote ASAS outputs:")
    print(f"  partial: {out['partial']} ({n} records)")
    print(f"  phased loci: {len(out['loci'])}")
    if args.ase:
        from .stages.asas import ase_sample
        if not args.ref_genome:
            sys.stderr.write("iflas asas --ase: needs --ref-genome and NGS reads (--ngs-left/--ngs-right).\n")
            return 2
        try:
            ase = ase_sample(args.base_dir, args.ref_genome, ngs_left=args.ngs_left, ngs_right=args.ngs_right,
                             paired=not args.ngs_single, ref_fa=args.ref_fa, alt_fa=args.alt_fa, use_freebayes=args.freebayes,
                             threads=args.threads or 1)
        except (FileNotFoundError, ValueError) as exc:
            sys.stderr.write(f"iflas asas --ase: {exc}\n")
            return 2
        print("Wrote ASE outputs:")
        for path in ase["merged"]:
            print(f"  {path} ({_count_nonempty_lines(path, missing_ok=True)} genes)")
    return 0


def _cmd_characterize(args) -> int:
    from .stages.characterize import characterize_sample

    try:
        out = characterize_sample(
            args.base_dir,
            out_dir=args.outdir,
            pa_grouped=args.pa_grouped,
            ngs=True if args.ngs else None,
            splice_site=not args.no_splice_site,
            ref_gpe=args.ref_gpe,
            ref_genome=args.ref_genome,
        )
    except (FileNotFoundError, ValueError) as exc:
        sys.stderr.write(f"iflas characterize: {exc}\n")
        return 2

    print("Wrote AS characterization outputs:")
    print(f"  output dir: {out['out_dir']}")
    for as_type, path in sorted(out["confident"].items()):
        n = _count_nonempty_lines(path, missing_ok=True)
        print(f"  confident.{as_type}.lst: {path} ({n} events)")
    if out["splice_site"]:
        print(f"  splice-site tables: {len(out['splice_site'])} AS types")
    for as_type, path in sorted(out.get("statistics", {}).items()):
        print(f"  statistic.{as_type}.tsv: {_count_nonempty_lines(path, missing_ok=True)} lines")
    return 0


def _cmd_go(args) -> int:
    from .stages.go import run_go_enrichment

    try:
        out = run_go_enrichment(
            args.target_gene_files,
            args.sample_names,
            args.gene2go_file,
            out_name=args.out_name,
            cutoff=args.cutoff,
            filter_by=args.filter_by,
            show_category=args.show_category,
            out_dir=args.outdir,
        )
    except (FileNotFoundError, ValueError) as exc:
        sys.stderr.write(f"iflas go: {exc}\n")
        return 2

    print("Ran GO enrichment via Rscript")
    print(f"  output dir: {out['go_dir']}")
    print(f"  table: {out['result_table']}")
    print(f"  plot: {out['result_pdf']}")
    return 0


def _cmd_diff_as(args) -> int:
    from .stages.diff_as import diff_as_project, run_rmats_from_bam_table

    try:
        if args.config:
            cfg = _load_config(args)
            out = diff_as_project(cfg, args.comp_cond, pu_filter=args.pu_filter, go=args.go, gene2go=args.gene2go,
                                  threads=args.threads)
        else:
            if not args.comp_cond or not args.gtf:
                sys.stderr.write("iflas diff_as: provide -c/-s (project mode) or --comp-cond and --gtf (table mode).\n")
                return 2
            out = run_rmats_from_bam_table(args.comp_cond, args.gtf, out_dir=args.outdir, conditions=args.conditions,
                                           threads=args.threads, sig=True)
    except (FileNotFoundError, ValueError) as exc:
        sys.stderr.write(f"iflas diff_as: {exc}\n")
        return 2

    print("Ran diff_as via rMATS-turbo")
    print(f"  output dir: {out['das_dir']}")
    for comp in out["comparisons"]:
        print(f"  {comp['name']}: {comp['out_dir']}" + (f" (significant: {Path(comp['sig']['SE']).parent})" if comp.get("sig") else ""))
    return 0


def _cmd_visual_as(args) -> int:
    from .stages.visual_as import visual_as_merged, visual_as_sample

    ngs = False if args.no_ngs else None
    try:
        if args.merge:
            outputs = visual_as_merged(args.base_dir, args.ref_gpe, genes=args.gene, gene_file=args.gene_file,
                                       out_dir=args.outdir, ngs=ngs, threads=args.threads or 1)
        else:
            if len(args.base_dir) != 1:
                sys.stderr.write("iflas visual_as: give one --base-dir, or several with --merge.\n")
                return 2
            outputs = visual_as_sample(args.base_dir[0], args.ref_gpe, genes=args.gene, gene_file=args.gene_file,
                                       out_dir=args.outdir, isoform_bed=args.isoform_bed, ngs=ngs,
                                       threads=args.threads or 1)
    except (FileNotFoundError, ValueError, KeyError) as exc:
        sys.stderr.write(f"iflas visual_as: {exc}\n")
        return 2

    print("Wrote visual_as outputs:")
    for out in outputs:
        print(f"  {out['gene']}: {out['pdf']} ({out['region']})")
    return 0


def _cmd_mapping(args) -> int:
    from . import orchestrator

    cfg = _load_config(args)
    orchestrator.restore_preproc_outputs(cfg)            # BUG-198: the results of an earlier `iflas preproc`
    if getattr(args, "junc_comb_support", None) is not None:
        cfg.minimap2.junc_comb_support = args.junc_comb_support
    selected_names = set(args.sample or [])
    if selected_names:
        samples = [s for s in cfg.samples if s.sample_name in selected_names]
        missing = sorted(selected_names - {s.sample_name for s in samples})
        if missing:
            sys.stderr.write("iflas mapping: unknown sample(s): " + ", ".join(missing) + "\n")
            return 2
    else:
        samples = list(cfg.samples)
    if not samples:
        sys.stderr.write("iflas mapping: config/sample sheet contains no samples.\n")
        return 2

    threads = max(1, int(args.threads or cfg.options.threads))
    results = []
    for sample in samples:
        job = orchestrator.SampleJob(
            stage="mapping",
            config=cfg,
            sample=sample,
            ref=cfg.ref_for(sample.ref_strain),
            base_dir=orchestrator.sample_base_dir(cfg, sample),
            threads=threads,
        )
        results.append((sample.sample_name, job.base_dir, orchestrator._stage_mapping(job)))

    print(f"Ran mapping for {len(results)} sample(s)")
    for sample_name, base_dir, outputs in results:
        print(f"  {sample_name}: {base_dir}/mapping")
        for key in sorted(outputs):
            print(f"    {key}: {outputs[key]}")
    return 0


def _cmd_preproc(args) -> int:
    from . import orchestrator

    cfg = _load_config(args)
    selected_names = set(args.sample or [])
    if selected_names:
        samples = [s for s in cfg.samples if s.sample_name in selected_names]
        missing = sorted(selected_names - {s.sample_name for s in samples})
        if missing:
            sys.stderr.write("iflas preproc: unknown sample(s): " + ", ".join(missing) + "\n")
            return 2
    else:
        samples = list(cfg.samples)
    if not samples:
        sys.stderr.write("iflas preproc: config/sample sheet contains no samples.\n")
        return 2

    threads = max(1, int(args.threads or cfg.options.threads))
    results = []
    for sample in samples:
        job = orchestrator.SampleJob(
            stage="preproc",
            config=cfg,
            sample=sample,
            ref=cfg.ref_for(sample.ref_strain),
            base_dir=orchestrator.sample_base_dir(cfg, sample),
            threads=threads,
        )
        results.append((sample.sample_name, job.base_dir, orchestrator._stage_preproc(job)))

    print(f"Ran preproc for {len(results)} sample(s)")
    for sample_name, base_dir, outputs in results:
        print(f"  {sample_name}: {base_dir}/preprocess")
        for key in sorted(outputs):
            if outputs[key]:
                print(f"    {key}: {outputs[key]}")
    return 0


def _cmd_collapse(args) -> int:
    from . import orchestrator

    cfg = _load_config(args)
    orchestrator.restore_preproc_outputs(cfg)            # BUG-198: the results of an earlier `iflas preproc`
    selected_names = set(args.sample or [])
    if selected_names:
        samples = [s for s in cfg.samples if s.sample_name in selected_names]
        missing = sorted(selected_names - {s.sample_name for s in samples})
        if missing:
            sys.stderr.write("iflas collapse: unknown sample(s): " + ", ".join(missing) + "\n")
            return 2
    else:
        samples = list(cfg.samples)
    if not samples:
        sys.stderr.write("iflas collapse: config/sample sheet contains no samples.\n")
        return 2

    threads = max(1, int(cfg.options.threads))
    results = []
    for sample in samples:
        reason = orchestrator.no_long_reads("collapse", sample)
        if reason:
            print(f"Skip collapse for {sample.sample_name}: {reason}")
            continue
        job = orchestrator.SampleJob(
            stage="collapse",
            config=cfg,
            sample=sample,
            ref=cfg.ref_for(sample.ref_strain),
            base_dir=orchestrator.sample_base_dir(cfg, sample),
            threads=threads,
        )
        results.append((sample.sample_name, job.base_dir, orchestrator._stage_collapse(job)))

    if not results:
        sys.stderr.write("iflas collapse: none of the selected samples has long reads.\n")
        return 2
    print(f"Ran collapse for {len(results)} sample(s)")
    for sample_name, base_dir, outputs in results:
        print(f"  {sample_name}: {base_dir}/collapse")
        for key in sorted(outputs):
            print(f"    {key}: {outputs[key]}")
    return 0


def _cmd_pu_filter(args) -> int:
    from .stages.pu_filter import filter_hq_isoforms
    out = filter_hq_isoforms(
        args.base_dir, args.ref_genome, args.ref_bed,
        sam_file=args.sam, junction_file=args.junctions,
        pos_fl_coverage=args.pos_fl_coverage, pos_min_junc_rpkm=args.pos_min_junc_rpkm,
        filter_score=args.filter_score, model=args.model, out_dir=args.outdir,
        feature_file=args.feature_file, random_state=args.random_state, select_best_model=args.select_best_model,
        auto_filter_score=args.auto_filter_score, draw_auc=args.draw_auc)
    n = _count_nonempty_lines(out["hq_bed"], missing_ok=True)
    print(f"Wrote {out['hq_bed']} ({n} high-quality isoforms; model {out['model']}, filter_score {out['filter_score']:.3f})")
    return 0


def _cmd_refine(args) -> int:
    from .stages.refine import refine_isoforms
    out = refine_isoforms(args.collapsed_gff, args.ref_genome, args.ref_gpe, args.ref_bed,
                          args.outdir or ".", ngs_junctions=args.ngs_junctions,
                          adjust=not args.no_strand_adjust)
    n = _count_nonempty_lines(out["isoformGrouped"])
    print(f"Wrote {out['isoformGrouped']} ({n} grouped isoforms)")
    print(f"  anno: {out['anno']}  novel: {out['novel']}")
    return 0


def _cmd_all(args) -> int:
    from . import orchestrator

    stages = tuple(s.strip() for s in args.stages.split(",") if s.strip()) if args.stages else None
    cfg = _load_config(args)
    results = orchestrator.run_pipeline(
        cfg,
        stages=stages,
        max_workers=args.jobs,
        merge_samples=cfg.options.merge_data_from_same_strain,
    )
    skipped = [r for r in results if "skipped" in r.outputs]
    print(f"Ran {len(results) - len(skipped)} stage/sample jobs" + (f", skipped {len(skipped)}" if skipped else ""))
    for result in results:
        note = f" (skipped: {result.outputs['skipped']})" if "skipped" in result.outputs else ""
        print(f"  {result.stage}: {result.sample_name} -> {result.base_dir}{note}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="iflas",
        description="iFLAS: integrated Full-Length Alternative Splicing analysis",
    )
    parser.add_argument("-V", "--version", action="version", version=f"iflas {__version__}")
    sub = parser.add_subparsers(title="command", dest="command", metavar="<command>")

    p_ref = sub.add_parser("init-ref",
                           help="Prepare the reference files and indexes from a genome FASTA and a GTF or GFF3.")
    p_ref.add_argument("-f", "--fasta", required=True, help="genome FASTA (.gz is decompressed into the output dir)")
    p_ref.add_argument("--gtf", help="gene annotation GTF (.gz accepted)")
    p_ref.add_argument("--gff", help="gene annotation GFF3, .gz accepted (converted via gffread if no --gtf)")
    p_ref.add_argument("-o", "--outdir", help="output dir for derived files (default: genome dir)")
    p_ref.add_argument("--hisat2", action="store_true", help="also build a HISAT2 index")
    p_ref.add_argument("--no-minimap2", action="store_true", help="skip the minimap2 index")
    p_ref.add_argument("-t", "--threads", type=int, default=1, help="threads for hisat2-build (default: 1)")
    p_ref.set_defaults(func=_cmd_init_ref)

    p_mig = sub.add_parser("migrate-config",
                           help="Convert an iFLAS 1 .cfg file into config.yaml and samples.tsv.")
    p_mig.add_argument("input", help="iFLAS 1 config file (.cfg)")
    p_mig.add_argument("-o", "--out", default="config.yaml", help="output YAML (default: config.yaml)")
    p_mig.add_argument("-s", "--samples", default="samples.tsv",
                       help="output sample-sheet TSV (default: samples.tsv)")
    p_mig.set_defaults(func=_cmd_migrate_config)

    from .orchestrator import DEFAULT_STAGES

    p_all = sub.add_parser("all",
                           help="Run the whole pipeline (preproc to palen_as) or selected stages.")
    p_all.add_argument("-c", "--config", required=True, help="config.yaml")
    p_all.add_argument("-s", "--samples", required=True, help="samples.tsv")
    p_all.add_argument("--stages", default=None,
                       help="comma-separated stages to run; default: the whole chain " + ",".join(DEFAULT_STAGES)
                            + ". In the default run a stage a sample lacks the inputs for is skipped with a notice "
                              "(ase without NGS, palen_as without a poly(A) table / nanopore fast5, diff_as without per-sample "
                              "alignments of two conditions); explicitly listed stages must run.")
    p_all.add_argument("-j", "--jobs", type=int, help="parallel sample jobs per stage (default: number of samples)")
    p_all.set_defaults(func=_cmd_all)

    p_pre = sub.add_parser("preproc",
                           help="Preprocess raw PacBio, Nanopore and NGS reads (ccs, lima, isoseq, fastp, fmlrc2).")
    p_pre.add_argument("-c", "--config", required=True, help="config.yaml")
    p_pre.add_argument("-s", "--samples", required=True, help="samples.tsv")
    p_pre.add_argument("--sample", action="append",
                       help="sample_name to run (repeatable; default: all samples)")
    p_pre.add_argument("-t", "--threads", type=int,
                       help="threads per sample (default: options.threads)")
    p_pre.set_defaults(func=_cmd_preproc)

    p_map = sub.add_parser("mapping",
                           help="Map long reads (minimap2) and NGS reads (HISAT2) to the genome; NGS junctions with regtools.")
    p_map.add_argument("-c", "--config", required=True, help="config.yaml")
    p_map.add_argument("-s", "--samples", required=True, help="samples.tsv")
    p_map.add_argument("--sample", action="append",
                       help="sample_name to run (repeatable; default: all samples)")
    p_map.add_argument("-t", "--threads", type=int,
                       help="threads per sample (default: options.threads)")
    p_map.add_argument("-jcs", "--junc-comb-support", dest="junc_comb_support", type=int,
                       help="FLNC reads sharing a junction chain needed to keep them "
                            "(default: minimap2.junc_comb_support = 2)")
    p_map.set_defaults(func=_cmd_mapping)

    p_col = sub.add_parser("collapse",
                           help="Collapse corrected reads into high-confidence isoforms (cDNA_Cupcake).")
    p_col.add_argument("-c", "--config", required=True, help="config.yaml")
    p_col.add_argument("-s", "--samples", required=True, help="samples.tsv")
    p_col.add_argument("--sample", action="append",
                       help="sample_name to run (repeatable; default: all samples)")
    p_col.set_defaults(func=_cmd_collapse)

    p_rfn = sub.add_parser("refine",
                           help="Correct isoform strands and junctions (with NGS junctions) and assign isoforms to genes.")
    p_rfn.add_argument("--collapsed-gff", required=True, help="cupcake collapsed GFF (tofu.collapsed.good.gff)")
    p_rfn.add_argument("--ref-genome", required=True, help="reference genome FASTA")
    p_rfn.add_argument("-r", "--ref-gpe", required=True, help="reference genePredExt (ref_gpe)")
    p_rfn.add_argument("--ref-bed", required=True, help="reference bed12+ with gene name (ref_bed)")
    p_rfn.add_argument("--ngs-junctions", help="NGS junctions.bed to enable junction refinement (optional)")
    p_rfn.add_argument("--no-strand-adjust", action="store_true", help="skip strand adjustment")
    p_rfn.add_argument("-o", "--outdir", help="output dir (writes refine/*; default: .)")
    p_rfn.set_defaults(func=_cmd_refine)

    p_puf = sub.add_parser("pu_filter",
                           help="Remove low-quality novel isoforms with PU learning.")
    p_puf.add_argument("-b", "--base-dir", required=True,
                       help="sample dir holding refine/ + collapse/ (+ mapping/ NGS) outputs")
    p_puf.add_argument("--ref-genome", required=True, help="reference genome FASTA")
    p_puf.add_argument("--ref-bed", required=True, help="reference bed12+ with gene name (ref_bed)")
    p_puf.add_argument("--sam", help="NGS reassembly tmp.bam (default: auto-detect under base-dir)")
    p_puf.add_argument("--junctions", help="NGS junctions.bed (default: auto-detect under base-dir)")
    p_puf.add_argument("-pos_fl_cov", "--pos-fl-coverage", dest="pos_fl_coverage", type=int, default=2,
                       help="min FL count for a positive (annotated) training isoform [2]")
    p_puf.add_argument("-pos_min_junc_rpkm", "--pos-min-junc-rpkm", dest="pos_min_junc_rpkm", type=float, default=0.05,
                       help="min junction RPKM for a positive training isoform [0.05]")
    p_puf.add_argument("-filter_score", "--filter-score", dest="filter_score", type=float, default=0.5,
                       help="PU score threshold to keep a novel isoform [0.5]")
    p_puf.add_argument("-select_best_model", "--select-best-model", dest="select_best_model", action="store_true", default=False,
                       help="pick the classifier with the best cross-validated/held-out AUC instead of GB")
    p_puf.add_argument("-auto_filter_score", "--auto-filter-score", dest="auto_filter_score", action="store_true", default=False,
                       help="with --select-best-model: threshold = 95th percentile of the negatives' PU scores")
    p_puf.add_argument("-draw_auc", "--draw-auc", dest="draw_auc", action="store_true", default=False,
                       help="with --select-best-model: write hqIsoforms/PU.bagging.ROC_AUC.pdf")
    p_puf.add_argument("--random-state", dest="random_state", type=int, default=0,
                       help="seed for the PU bagging and the training-set shuffle (reproducible scores)")
    p_puf.add_argument("--model", default="GB", choices=["RF", "GB", "DT", "SVM", "NB"],
                       help="base classifier for the PU bagging ensemble [GB]")
    p_puf.add_argument("-feature_file", "--feature-file", dest="feature_file",
                       help="precomputed isoFeatures.txt table to reuse instead of recalculating features")
    p_puf.add_argument("-o", "--outdir", help="output dir (writes hqIsoforms/*; default: base-dir)")
    p_puf.set_defaults(func=_cmd_pu_filter)

    p_fas = sub.add_parser("find_as",
                           help="Identify AS events (IR, SE, A3SS, A5SS) from the high-confidence isoforms.")
    p_fas.add_argument("-b", "--base-dir",
                       help="sample dir with hqIsoforms/hq.collapsed.bed12+ and collapse/tofu.collapsed.group.txt")
    p_fas.add_argument("--isoform-bed",
                       help="isoform bed12+ (col15 = group key); default: <base-dir>/hqIsoforms/hq.collapsed.bed12+")
    p_fas.add_argument("-g", "--group",
                       help="cupcake tofu.collapsed.group.txt (default: <base-dir>/collapse/tofu.collapsed.group.txt)")
    p_fas.add_argument("-r", "--ref-gpe", required=True,
                       help="reference genePredExt (ref_gpe) for annotation")
    p_fas.add_argument("-o", "--outdir", help="output dir (writes LR/*.bed; default: .)")
    p_fas.add_argument("--ngs-junctions",
                       help="NGS junctions.bed for confident filtering (base-dir mode auto-detects mapping/rna-seq/reassembly/junctions.bed)")
    p_fas.set_defaults(func=_cmd_find_as)

    p_fpa = sub.add_parser("find_pa",
                           help="Identify poly(A) sites and APA in a sample (-b), or only cluster poly(A) sites.")
    p_fpa.add_argument("-b", "--base-dir",
                       help="sample dir with refine/reads.assigned.unambi.bed12+, refine/isoformGrouped.bed12+ and "
                            "collapse/tofu.collapsed.group.txt; writes as_events/pa/ (PA.bed6+, paGrouped tables, motif/)")
    p_fpa.add_argument("--ref-genome", help="reference genome FASTA (base-dir mode)")
    p_fpa.add_argument("--ref-size", help="chromosome sizes (ref_size / genome.len.txt; base-dir mode)")
    p_fpa.add_argument("-t", "--threads", type=int, help="threads (base-dir mode)")
    p_fpa.add_argument("--reads-bed",
                       help="reads.assigned.unambi.bed12+ from refine (clustering-only mode)")
    p_fpa.add_argument("-g", "--group",
                       help="tofu.collapsed.group.txt from collapse (clustering-only mode)")
    p_fpa.add_argument("-o", "--out",
                       help="output paCluster.bed8+ path (clustering-only mode)")
    p_fpa.add_argument("-pa_sup", "--pa-support", dest="pa_support", type=int, default=5,
                       help="minimum FL read support for a PA cluster [5]")
    p_fpa.add_argument("-pa_rpkm", "--pa-rpkm", dest="pa_rpkm", type=float, default=0,
                       help="RPKM(PAC) threshold [0]")
    p_fpa.add_argument("--pa-dist", type=int, default=24,
                       help="distance used to merge nearby PAC peaks [24]")
    p_fpa.add_argument("--window-size", type=int, default=3,
                       help="window size for PAC peak detection [3]")
    p_fpa.add_argument("-conf_pa", "--confident-pa", dest="confident_pa",
                       help="optional confident PA/PAC file for Known/Novel annotation")
    p_fpa.set_defaults(func=_cmd_find_pa)

    p_char = sub.add_parser("characterize",
                            help="Summarize AS events: annotated and confident event lists, splice-site tables and "
                                 "comparison with the reference.")
    p_char.add_argument("-b", "--base-dir", required=True,
                        help="sample dir with as_events/ordinary_as/ (find_as output)")
    p_char.add_argument("-o", "--outdir",
                        help="output dir (default: <base-dir>/as_events/characterization)")
    p_char.add_argument("--pa-grouped",
                        help="reads.paGrouped.tsv for the APA lists (default: <base-dir>/as_events/pa/reads.paGrouped.tsv when present)")
    p_char.add_argument("--ngs", action="store_true",
                        help="force the NGS branch (Common/LR.specific lists + Venn plots; default: auto-detect ordinary_as/NGS)")
    p_char.add_argument("-r", "--ref-gpe",
                        help="reference genePredExt: enables the back half (reference AS classification, anno/novel "
                             "statistics, APA comparison, distance to TTS)")
    p_char.add_argument("--ref-genome", help="reference genome FASTA for the splice-site dinucleotide profiles")
    p_char.add_argument("--no-splice-site", action="store_true",
                        help="only assemble the event lists; skip the splice-site decompositions")
    p_char.set_defaults(func=_cmd_characterize)

    p_vas = sub.add_parser("visual_as",
                           help="Plot gene structures: reference models, iFLAS isoforms, NGS coverage and junctions (pyGenomeTracks).")
    p_vas.add_argument("-b", "--base-dir", required=True, action="append",
                       help="sample dir with hqIsoforms/hq.collapsed.bed12+ (repeat with --merge to stack samples)")
    p_vas.add_argument("--merge", action="store_true",
                       help="one plot per gene with every given sample stacked (written to isoViewer_sample_merged/)")
    p_vas.add_argument("--no-ngs", action="store_true", help="skip the NGS coverage / junction tracks")
    p_vas.add_argument("-t", "--threads", type=int, help="genes plotted in parallel (per-sample mode)")
    p_vas.add_argument("-r", "--ref-gpe", required=True,
                       help="reference genePredExt (ref_gpe) for reference gene tracks")
    p_vas.add_argument("-g", "--gene", action="append",
                       help="target gene; repeat or comma-separate. Default: all genes in the isoform BED")
    p_vas.add_argument("--gene-file",
                       help="target gene file, one gene id per line")
    p_vas.add_argument("--isoform-bed",
                       help="hq isoform bed12+ (default: <base-dir>/hqIsoforms/hq.collapsed.bed12+)")
    p_vas.add_argument("-o", "--outdir",
                       help="output root (default: <base-dir>/visual_as)")
    p_vas.set_defaults(func=_cmd_visual_as)

    p_diff = sub.add_parser("diff_as",
                            help="Differential AS between conditions with rMATS (from config.yaml and samples.tsv, or a BAM table).")
    p_diff.add_argument("-c", "--config", help="config.yaml (project mode: merges the samples' AS isoforms and uses their NGS BAMs)")
    p_diff.add_argument("-s", "--samples", help="samples.tsv (project mode)")
    p_diff.add_argument("-d", "--comp-cond", dest="comp_cond",
                        help="7-column BAM table (project/sampleName/condition/repeat/bamFile/paired/readsLength) or, in "
                             "project mode, a comma-separated list of conditions to compare (default: all)")
    p_diff.add_argument("--pu-filter", dest="pu_filter", action="store_true",
                        help="merge the PU-filtered HQ isoforms instead of the isoforms with confident AS events (-pu_filter in iFLAS 1)")
    p_diff.add_argument("--go", action="store_true", help="GO enrichment of the differentially spliced genes (needs gene2go)")
    p_diff.add_argument("--gene2go", help="gene -> GO mapping file (default: options.gene2go)")
    p_diff.add_argument("--gtf",
                        help="merged isoform GTF used by rMATS (table mode)")
    p_diff.add_argument("-o", "--outdir", default=".",
                        help="base output dir; writes <outdir>/das [.]")
    p_diff.add_argument("--conditions",
                        help="comma-separated condition subset/order; default: all table conditions")
    p_diff.add_argument("-t", "--threads", type=int, default=1,
                        help="threads for --tstat and --nthread [1]")
    p_diff.set_defaults(func=_cmd_diff_as)

    p_asas = sub.add_parser("asas",
                            help="Identify allele-specific AS: phase the full-length reads of each locus with IsoPhase "
                                 "(--ref-genome and --flnc) and test the haplotypes against AS isoform pairs.")
    p_asas.add_argument("-b", "--base-dir", required=True,
                        help="sample dir with hqIsoforms/, collapse/, as_events/ (and asas/by_loci/ when not phasing here)")
    p_asas.add_argument("--ref-genome", help="reference genome FASTA: run IsoPhase (select_loci_to_phase + run_phaser) first")
    p_asas.add_argument("--flnc", help="FLNC reads (fasta/fastq) for IsoPhase")
    p_asas.add_argument("--collapsed-gff", help="cupcake collapsed GFF (default: <base-dir>/collapse/tofu.collapsed.good.gff)")
    p_asas.add_argument("--min-coverage", type=int, default=10, help="minimum FLNC coverage of a locus to phase [10]")
    p_asas.add_argument("-t", "--threads", type=int, help="loci phased in parallel")
    p_asas.add_argument("--ase", action="store_true",
                        help="also quantify allele-specific expression of the phased loci from the NGS reads (-ase in iFLAS 1)")
    p_asas.add_argument("--ngs-left", help="NGS left reads (files of one replicate joined by ',', replicates separated by ';')")
    p_asas.add_argument("--ngs-right", help="NGS right reads (same syntax)")
    p_asas.add_argument("--ngs-single", action="store_true", help="NGS reads are single-end")
    p_asas.add_argument("--no-ngs-confirm", action="store_true",
                        help="keep phased loci without an NGS-confirmed SNP (as in iFLAS 1)")
    p_asas.add_argument("--freebayes", action="store_true", help="confirm heterozygous SNPs with freebayes on the NGS BAM (-fbs in iFLAS 1)")
    p_asas.add_argument("--ref-fa", help="pre-built REF haplotype transcript FASTA for salmon (-ref_fa in iFLAS 1)")
    p_asas.add_argument("--alt-fa", help="pre-built ALT haplotype transcript FASTA for salmon (-alt_fa in iFLAS 1)")
    p_asas.add_argument("-o", "--outdir",
                        help="output dir (default: <base-dir>/asas)")
    p_asas.add_argument("--isoform-bed",
                        help="hq isoform bed12+ (default: <base-dir>/hqIsoforms/hq.collapsed.bed12+)")
    p_asas.add_argument("-g", "--group",
                        help="tofu.collapsed.group.txt (default: <base-dir>/collapse/tofu.collapsed.group.txt)")
    p_asas.add_argument("--as-dir",
                        help="directory of the confident long-read AS events (default: <base-dir>/as_events/ordinary_as/LR)")
    p_asas.add_argument("--loci-dir",
                        help="IsoPhase by_loci dir (default: <outdir>/by_loci)")
    p_asas.add_argument("--filter-by-count", type=int, default=2,
                        help="minimum FL read support per isoform [2]")
    p_asas.add_argument("--merge-by-junc", action="store_true",
                        help="merge AS isoforms by junction-combination before testing")
    p_asas.set_defaults(func=_cmd_asas)

    p_pal = sub.add_parser("palen_as",
                           help="Compare poly(A) tail lengths between AS isoforms and between APA sites (nanopolish or dorado tails).")
    p_pal.add_argument("-b", "--base-dir", required=True,
                       help="sample dir with hqIsoforms/, collapse/, as_events/, and mapping/ outputs")
    p_pal.add_argument("--palen", action="append",
                       help="nanopolish polya_results.tsv file; repeat or comma-separate for multiple files")
    p_pal.add_argument("--fast5-dir", help="ONT fast5 directory: run nanopolish index/polya here (needs --ref-genome and preprocess/nanopore/rawFlnc.fq)")
    p_pal.add_argument("--ref-genome", help="reference genome FASTA for nanopolish polya")
    p_pal.add_argument("--sequencing-summary", help="basecalling sequencing_summary.txt for nanopolish index")
    p_pal.add_argument("-t", "--threads", type=int, help="nanopolish threads")
    p_pal.add_argument("-o", "--outdir",
                       help="output dir (default: <base-dir>/palenAS)")
    p_pal.add_argument("--isoform-bed",
                       help="hq isoform bed12+ (default: <base-dir>/hqIsoforms/hq.collapsed.bed12+)")
    p_pal.add_argument("-g", "--group",
                       help="tofu.collapsed.group.txt (default: <base-dir>/collapse/tofu.collapsed.group.txt)")
    p_pal.add_argument("--reads-bed",
                       help="assigned FLNC reads bed12+ (default: <base-dir>/mapping/flnc.addCVandID.bed12+)")
    p_pal.add_argument("--as-dir",
                       help="directory of the confident long-read AS events (default: <base-dir>/as_events/ordinary_as/LR)")
    p_pal.add_argument("--pa-cluster",
                       help="paCluster.bed8+ for APA tail-length testing (default: <base-dir>/as_events/pa/paCluster.bed8+)")
    p_pal.add_argument("-pa_sup", "--pa-support", dest="pa_support", type=int, default=10,
                       help="minimum FL read support for AS/APA tail-length tests [10]")
    p_pal.add_argument("-conf_pac", "--confident-pac", dest="confident_pac", action="append",
                       help="confident PAC file(s) used to filter APA results; repeat or comma-separate")
    p_pal.add_argument("--ref-bed",
                       help="reference bed12+ required when --confident-pac is used")
    p_pal.add_argument("--no-merge-by-junc", action="store_true",
                       help="do not merge AS isoforms by junction-combination before testing")
    p_pal.set_defaults(func=_cmd_palen_as)

    p_go = sub.add_parser("go",
                          help="GO enrichment analysis and plots for one or more gene sets.")
    p_go.add_argument("-tg", "--target-genes", dest="target_gene_files", required=True,
                      help="target gene file, or comma-separated files for multiple samples")
    p_go.add_argument("-s", "--sample-name", dest="sample_names", required=True,
                      help="sample name, or comma-separated names matching --target-genes")
    p_go.add_argument("-bg", "--gene2go", dest="gene2go_file", required=True,
                      help="gene-to-GO mapping TSV used for enrichment")
    p_go.add_argument("-o", "--out-name", dest="out_name", default="goEnrichment",
                      help="output prefix [goEnrichment]")
    p_go.add_argument("--outdir", default=".", help="base output dir; writes <outdir>/GO [.]")
    p_go.add_argument("-cutoff", "--cutoff", type=float, default=0.05,
                      help="cutoff used to filter enrichment results [0.05]")
    p_go.add_argument("-filterBy", "--filter-by", dest="filter_by",
                      choices=["pvalue", "p.adjust"], default="p.adjust",
                      help="column used for filtering [p.adjust]")
    p_go.add_argument("-showCategory", "--show-category", dest="show_category", type=int, default=20,
                      help="number of GO categories shown in the plot [20]")
    p_go.set_defaults(func=_cmd_go)

    p_report = sub.add_parser("report",
                              help="Make plots and an HTML report from the finished stages.")
    p_report.add_argument("-c", "--config",
                          help="config.yaml used to locate output directories and samples")
    p_report.add_argument("-s", "--samples",
                          help="samples.tsv used to locate output directories and samples")
    p_report.add_argument("--sample", action="append",
                          help="sample_name to include (repeatable; default: all samples)")
    p_report.add_argument("-all", dest="all", action="store_true", default=False,
                          help="all plots and the HTML report")
    p_report.add_argument("-basic", dest="basic", action="store_true", default=False,
                          help="basic read and isoform statistics")
    p_report.add_argument("-asp", dest="asp", action="store_true", default=False,
                          help="AS characterization plots")
    p_report.add_argument("-geneStruc", dest="gene_struc", action="store_true", default=False,
                          help="merge the visual_as gene-structure PDFs")
    p_report.add_argument("-asas", dest="asas", action="store_true", default=False,
                          help="allele-specific AS plots")
    p_report.add_argument("-palen", dest="palen", action="store_true", default=False,
                          help="poly(A) tail length plots for AS isoforms and APA sites")
    p_report.add_argument("-diff", dest="diff", action="store_true", default=False,
                          help="differential AS plots (das/*.sigDiffAS)")
    p_report.add_argument("-html", "--html", dest="html", action="store_true", default=False,
                          help="write reports/html and reports/json/search.json")
    p_report.add_argument("--novel-hq-as", dest="novel_hq_as", action="store_true", default=False,
                          help="score plots of the novel HQ isoforms (isoIdentity/isoformScore.txt)")
    p_report.set_defaults(func=_cmd_report)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "command", None):
        parser.print_help()
        return 1
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
