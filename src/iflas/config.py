"""Configuration system: YAML (global params) + TSV (sample sheet), validated by pydantic.

Replaces the legacy INI (`Config.py`):
  * one `config.yaml` for reference + tool params + runtime,
  * one `samples.tsv` (one row per sample/replicate — no more `;`/`,` ambiguity),
  * pydantic validation with clear errors and defaults,
  * `${ENV}` expansion and relative-path resolution (relative to the config file's dir).

Field names mirror the legacy tags so behavior carries over; `migrate_ini` converts old configs.
"""

from __future__ import annotations

import configparser
import csv
import os
from typing import Callable, Optional

import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator

# ---- file-valued fields (resolved to absolute paths during load) ----
_REF_FILE_FIELDS = {"ref_genome", "ref_size", "ref_gtf", "ref_gff", "ref_gpe", "ref_bed",
                    "ref_mm2_index", "hisat2_index"}
_SAMPLE_FILE_FIELDS = {"data_location", "primer", "data_processed_location", "polya_location", "fast5_dir",
                       "ngs_left_reads", "ngs_right_reads", "ngs_junctions"}
_NONE_TOKENS = {"", "none", "null", "na", "."}


def resolve_path(value: Optional[str], base_dir: str) -> Optional[str]:
    """Expand ~ and ${ENV}, then make relative paths absolute against base_dir."""
    if value is None:
        return None
    value = os.path.expandvars(os.path.expanduser(str(value)))
    if value and not os.path.isabs(value):
        value = os.path.normpath(os.path.join(base_dir, value))
    return value


class _Base(BaseModel):
    model_config = ConfigDict(extra="forbid")


class RefSection(_Base):
    ref_strain: str
    ref_genome: str
    ref_gtf: str

    @model_validator(mode="before")
    @classmethod
    def _gff_needs_init_ref(cls, data):
        # legacy converted ref_gff at start-up; here the GTF and the files derived from it come from `iflas init-ref`
        if isinstance(data, dict) and not data.get("ref_gtf") and data.get("ref_gff"):
            raise ValueError("ref_gtf is required: convert the GFF3 and derive the reference files with "
                             f"`iflas init-ref -f <genome.fa> --gff {data['ref_gff']} -o <dir>`, then use the printed ref_* paths")
        return data
    ref_size: Optional[str] = None
    ref_gff: Optional[str] = None
    ref_gpe: Optional[str] = None
    ref_bed: Optional[str] = None
    ref_mm2_index: Optional[str] = None
    hisat2_index: Optional[str] = None

    @model_validator(mode="after")
    def _uncompressed(self):
        # samtools faidx, hisat2-build, freebayes and the annotation readers need plain files (BUG-194)
        from .refprep import is_gzipped
        for field in ("ref_genome", "ref_gtf"):
            path = getattr(self, field)
            if os.path.isfile(path) and is_gzipped(path):
                raise ValueError(f"{field} is compressed ({path}): `iflas init-ref` writes uncompressed copies of the "
                                 "genome and annotation; use the ref_* paths it prints")
        return self


class CcsParams(_Base):
    min_read_length: int = 50
    min_read_score: float = 0.75
    min_subread_length: int = 50
    min_ccs_length: int = 50
    min_predicted_accuracy: float = 0.9
    min_pass: int = 2


class Minimap2Params(_Base):
    mm2_index: Optional[str] = None
    max_intron_length: int = 50000
    junc_comb_support: int = 2   # legacy `mapping -jcs`: FLNC reads sharing a junction chain needed to keep them


class CollapseParams(_Base):
    min_identity: float = 0.9
    min_coverage: float = 0.9
    fl_coverage: int = 2
    max_fuzzy_junction: int = 5
    max_5_diff: int = 1000
    max_3_diff: int = 500
    dun_merge_5_shorter: bool = True


class PuFilterParams(_Base):
    filter_score: float = 0.5
    random_state: int = 0
    draw_auc: bool = False
    pos_fl_coverage: int = 2
    pos_min_junc_rpkm: float = 0.05
    select_best_model: bool = False
    auto_filter_score: bool = False
    feature_file: Optional[str] = None


class Options(_Base):
    threads: int = 4
    memory: str = "4000M"
    gene2go: Optional[str] = None
    merge_data_from_same_strain: bool = True
    # cut the basecalled 3' poly(A) of ONT direct-RNA reads before mapping (BUG-178, adopted at CP-ONT 2026-10-15);
    # false reproduces the published pipeline
    ont_trim_polya: bool = True
    # allele-specific AS/ASE only on phased loci with an SNP the NGS reads confirm (BUG-201, adopted 2026-10-16);
    # false reproduces the published pipeline
    asas_ngs_confirm: bool = True


class Dirs(_Base):
    out_dir: str = "outDir"
    tmp_dir: Optional[str] = None


class Sample(_Base):
    sample_name: str
    ref_strain: str
    project_name: str = "iflas"
    strain: Optional[str] = None
    condition: Optional[str] = None
    replicate: Optional[str] = None
    tgs_plat: Optional[str] = None
    strategy: Optional[str] = None
    data_location: Optional[str] = None
    primer: Optional[str] = None
    data_processed_location: Optional[str] = None
    polya_location: Optional[str] = None
    fast5_dir: Optional[str] = None          # raw ONT signal directory for nanopolish when data_location holds basecalled reads
    ngs_left_reads: Optional[str] = None
    ngs_right_reads: Optional[str] = None
    ngs_reads_paired: str = "paired"
    ngs_reads_length: Optional[int] = None
    ngs_junctions: Optional[str] = None
    ngs_library_strand: Optional[int] = None
    use_fmlrc2: bool = True
    single_run_threads: int = 4
    flowcell_type: Optional[str] = None
    kit_type: Optional[str] = None


class MergedTgsSample(_Base):
    """Config-layer port of legacy ``commonObjs.MergedTgsSample``."""

    project_name: Optional[str] = None
    sample_name: Optional[str] = None
    ref_strain: Optional[str] = None
    strain: Optional[str] = None
    condition: Optional[str] = None

    tgs_plat: Optional[str] = None
    strategy: Optional[str] = None
    data_location: Optional[str] = None
    data_processed_location: list[Optional[str]] = Field(default_factory=list)
    polya_location: Optional[dict[str, str]] = None
    fast5_dir: Optional[dict[str, str]] = None

    ngs_left_reads: Optional[str] = None
    ngs_right_reads: Optional[str] = None
    ngs_reads_paired: str = "paired"
    ngs_reads_length: Optional[int] = None
    ngs_junctions: Optional[str] = None

    use_fmlrc2: bool = True
    single_run_threads: int = 1

    def __str__(self) -> str:
        return "%s:%s:%s, %s:%s:%s, %s:%s:%s:%s:%s, %d" % (
            self.project_name,
            self.sample_name,
            self.ref_strain,
            self.tgs_plat,
            self.strategy,
            self.data_processed_location,
            self.ngs_left_reads,
            self.ngs_right_reads,
            self.ngs_reads_paired,
            self.ngs_reads_length,
            self.ngs_junctions,
            self.single_run_threads,
        )

    __repr__ = __str__


class IflasConfig(_Base):
    reference: list[RefSection]
    ccs: CcsParams = CcsParams()
    minimap2: Minimap2Params = Minimap2Params()
    collapse: CollapseParams = CollapseParams()
    pu_filter: PuFilterParams = PuFilterParams()
    options: Options = Options()
    dirs: Dirs = Dirs()
    samples: list[Sample] = []

    def ref_for(self, strain: str) -> RefSection:
        for r in self.reference:
            if r.ref_strain == strain:
                return r
        raise KeyError(f"no reference for ref_strain {strain!r}")

    @model_validator(mode="after")
    def _check_ref_strains(self):
        known = {r.ref_strain for r in self.reference}
        for s in self.samples:
            if s.ref_strain not in known:
                raise ValueError(
                    f"sample {s.sample_name!r} references unknown ref_strain {s.ref_strain!r}; "
                    f"known: {sorted(known)}"
                )
        return self


def _merge_repeat_reads(samples: list[Sample], attr: str) -> Optional[str]:
    """Repeat i of the merged sample = repeat i of every row that has NGS reads, comma-joined (legacy ``mergeSample``).
    Rows without NGS reads contribute none (legacy crashed on them: BUG-176, e.g. an ONT row in a strain whose
    PacBio row has NGS)."""
    present = [reads for reads in (getattr(sample, attr) for sample in samples) if reads is not None]
    if not present:
        return None
    min_repeats = min(len(reads.split(";")) for reads in present)
    return ";".join(",".join(reads.split(";")[i] for reads in present) for i in range(min_repeats))


def _merge_strain(project_name, ref_strain, strain, data_objs, validate_file) -> MergedTgsSample:
    """The ``<strain>_all`` sample of one project / reference strain / strain group."""
    merged = MergedTgsSample()
    merged.project_name = project_name
    merged.sample_name = "{}_all".format(strain)
    merged.ref_strain = ref_strain
    merged.strain = strain
    merged.condition = "{}_all".format(strain)
    long_read = [sample for sample in data_objs if sample.tgs_plat or sample.data_location or sample.data_processed_location]
    platform = next((sample for sample in long_read if sample.tgs_plat), data_objs[0])
    merged.tgs_plat = platform.tgs_plat               # legacy took the first row's: None when an NGS-only row came first
    merged.strategy = platform.strategy
    merged.data_processed_location = [sample.data_processed_location for sample in long_read]   # NGS-only rows have none
    merged.polya_location = {sample.sample_name: sample.polya_location for sample in data_objs
                             if sample.polya_location is not None and validate_file(sample.polya_location)} or None
    merged.fast5_dir = {sample.sample_name: sample.fast5_dir for sample in data_objs
                        if sample.fast5_dir and os.path.isdir(sample.fast5_dir)} or None
    merged.ngs_left_reads = _merge_repeat_reads(data_objs, "ngs_left_reads")
    merged.ngs_right_reads = _merge_repeat_reads(data_objs, "ngs_right_reads")
    merged.ngs_reads_paired = "single" if bool(merged.ngs_left_reads) != bool(merged.ngs_right_reads) else "paired"
    merged.use_fmlrc2 = data_objs[0].use_fmlrc2
    return merged


def merge_samples_by_strain(
    samples: list[Sample],
    *,
    validate_file: Callable[[str], bool] | None = None,
) -> dict[str, dict[str, dict[Optional[str], MergedTgsSample]]]:
    """Merge samples by project/ref_strain/strain using legacy ``mergeSample`` rules."""

    validate_file = validate_file or os.path.isfile
    grouped: dict[str, dict[str, dict[Optional[str], list[Sample]]]] = {}
    for sample in samples:
        grouped.setdefault(sample.project_name, {}).setdefault(sample.ref_strain, {}).setdefault(
            sample.strain, []
        ).append(sample)

    merged_to_process: dict[str, dict[str, dict[Optional[str], MergedTgsSample]]] = {}
    for project_name, by_ref in grouped.items():
        for ref_strain, by_strain in by_ref.items():
            for strain, data_objs in by_strain.items():
                merged_to_process.setdefault(project_name, {}).setdefault(ref_strain, {})[strain] = _merge_strain(
                    project_name, ref_strain, strain, data_objs, validate_file)
    return merged_to_process


def _clean(value: Optional[str]) -> Optional[str]:
    "map empty/none-ish tokens to None"
    if value is None:
        return None
    return None if str(value).strip().lower() in _NONE_TOKENS else value


def read_samples_tsv(path: str, base_dir: Optional[str] = None) -> list[dict]:
    """Read a TSV sample sheet into a list of dicts (none-ish -> None; paths resolved)."""
    base_dir = base_dir if base_dir is not None else os.path.dirname(os.path.abspath(path))
    rows: list[dict] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for raw in reader:
            row = {k.strip(): _clean(v) for k, v in raw.items() if k and k.strip()}
            for f in _SAMPLE_FILE_FIELDS:
                if row.get(f):
                    row[f] = resolve_path(row[f], base_dir)
            rows.append(row)
    return rows


def load_config(config_path: str, samples_path: Optional[str] = None) -> IflasConfig:
    """Load config.yaml (+ optional samples.tsv) into a validated IflasConfig."""
    base_dir = os.path.dirname(os.path.abspath(config_path))
    with open(config_path) as fh:
        data = yaml.safe_load(fh) or {}

    # reference may be a single mapping or a list of mappings
    ref = data.get("reference", [])
    if isinstance(ref, dict):
        ref = [ref]
    for r in ref:
        for f in _REF_FILE_FIELDS:
            if r.get(f):
                r[f] = resolve_path(r[f], base_dir)
    data["reference"] = ref

    for key in ("out_dir", "tmp_dir"):     # BUG-199: were taken relative to the working directory, not this file
        if (data.get("dirs") or {}).get(key):
            data["dirs"][key] = resolve_path(data["dirs"][key], base_dir)
    if data.get("options", {}).get("gene2go"):
        data["options"]["gene2go"] = resolve_path(data["options"]["gene2go"], base_dir)
    if data.get("pu_filter", {}).get("feature_file"):
        data["pu_filter"]["feature_file"] = resolve_path(data["pu_filter"]["feature_file"], base_dir)

    # samples: from TSV if given, else inline 'samples' in the YAML
    if samples_path:
        data["samples"] = read_samples_tsv(samples_path)
    elif data.get("samples") is None:
        data["samples"] = []

    return IflasConfig.model_validate(data)


# ---------------------------------------------------------------------------
# Migration: legacy INI (Config.py) -> config.yaml + samples.tsv
# ---------------------------------------------------------------------------
_INI_SECTION_TO_BLOCK = {
    "ccscalling": "ccs",
    "minimap2": "minimap2",
    "collapse": "collapse",
    "pu_filter": "pu_filter",
    "hqisofilter": "pu_filter",  # legacy example uses section_type=hqIsoFilter (Config.py canonical: pu_filter)
    "optiontools": "options",
    "dirs": "dirs",
}


def migrate_ini(ini_path: str):
    """Parse a legacy INI config into (config_dict, samples_list, skipped_sections).

    Routes sections by their ``section_type`` value (as the legacy parser did).
    `buildInTools` (bash_location, ...) is intentionally dropped.
    """
    parser = configparser.ConfigParser(interpolation=None)  # tolerate % in paths
    parser.read(ini_path)
    config: dict = {"reference": []}
    samples: list[dict] = []
    skipped: list[str] = []
    for sec in parser.sections():
        items = dict(parser.items(sec))
        stype = items.pop("section_type", "").strip().lower()
        if stype == "refsection":
            config["reference"].append(items)
        elif stype == "datacfg":
            items.setdefault("sample_name", sec)
            samples.append(items)
        elif stype in _INI_SECTION_TO_BLOCK:
            config[_INI_SECTION_TO_BLOCK[stype]] = items
        else:  # buildInTools / unknown -> drop
            skipped.append(sec)
    return config, samples, skipped


def write_samples_tsv(samples: list[dict], path: str) -> None:
    """Write sample dicts to a TSV (sample_name, ref_strain first; missing -> 'None')."""
    cols: list[str] = []
    for s in samples:
        for k in s:
            if k not in cols:
                cols.append(k)
    front = [c for c in ("sample_name", "ref_strain") if c in cols]
    ordered = front + [c for c in cols if c not in front]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=ordered, delimiter="\t", restval="None")
        writer.writeheader()
        for s in samples:
            writer.writerow(s)


def migrate_config_file(ini_path: str, yaml_out: str, samples_out: str):
    """Convert a legacy INI into config.yaml + samples.tsv. Returns (config, samples, skipped)."""
    config, samples, skipped = migrate_ini(ini_path)
    with open(yaml_out, "w") as fh:
        yaml.safe_dump(config, fh, sort_keys=False, allow_unicode=True)
    write_samples_tsv(samples, samples_out)
    return config, samples, skipped
