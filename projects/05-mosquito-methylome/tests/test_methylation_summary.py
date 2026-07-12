from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]


def load_summary():
    path = PROJECT / "scripts/07_summarize_methylation_results.py"
    spec = importlib.util.spec_from_file_location("methylation_summary", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_bedmethyl_parser_and_malformed_count(tmp_path: Path) -> None:
    summary = load_summary()
    bed = tmp_path / "sample.bedmethyl"
    bed.write_text(
        "mtDNA\t9\t10\ta\t0\t+\t9\t10\t0,0,0\t25\t40.0\t10\t15\t0\t0\t0\t0\t0\n"
        "malformed\trow\n",
        encoding="utf-8",
    )
    sites, malformed = summary.read_sites(bed)
    assert malformed == 1
    assert len(sites) == 1
    assert sites[0].position_1based == 10
    assert sites[0].coverage == 25
    assert sites[0].fraction_modified == 0.4


def test_empty_motif_means_not_detected(tmp_path: Path) -> None:
    summary = load_summary()
    motif = tmp_path / "motif.tsv"
    motif.write_text("motif\tscore\n", encoding="utf-8")
    assert summary.motif_summary(motif) == ("not_detected", 0, "")
