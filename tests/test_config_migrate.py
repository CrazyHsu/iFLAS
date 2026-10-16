"""Migration test: legacy INI (test_data/config.cfg) -> config.yaml + samples.tsv, then reload."""

import os

from iflas.config import load_config, migrate_config_file

LEGACY_INI = os.path.join(os.path.dirname(__file__), "..", "test_data", "config.cfg")


def test_migrate_legacy_config(tmp_path):
    yaml_out = tmp_path / "config.yaml"
    tsv_out = tmp_path / "samples.tsv"
    config, samples, skipped = migrate_config_file(LEGACY_INI, str(yaml_out), str(tsv_out))

    # reference + samples extracted
    assert config["reference"] and config["reference"][0]["ref_strain"] == "B73"
    assert any(s["sample_name"] == "B73" for s in samples)
    # tool-param blocks captured (incl. hqIsoFilter -> pu_filter)
    assert "ccs" in config and "pu_filter" in config and "dirs" in config
    # buildInTools dropped
    assert "buildInTools" not in config

    # the migrated config reloads & validates
    cfg = load_config(str(yaml_out), str(tsv_out))
    known = {r.ref_strain for r in cfg.reference}
    assert cfg.samples and all(s.ref_strain in known for s in cfg.samples)
