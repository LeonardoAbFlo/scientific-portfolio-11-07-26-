from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]


def load_script(name: str):
    path = PROJECT / "scripts" / name
    spec = importlib.util.spec_from_file_location(name.replace(".", "_"), path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_genus_parses_gtdb_and_species_names() -> None:
    filtering = load_script("06_filter_contigs.py")
    assert filtering.genus("d__Bacteria;p__Pseudomonadota;g__Klebsiella;s__Klebsiella pneumoniae") == "Klebsiella"
    assert filtering.genus("Pseudomonas aeruginosa") == "Pseudomonas"
    assert filtering.genus("Unknown") == ""


def test_split_and_merge_round_trip(tmp_path: Path) -> None:
    assemblies = tmp_path / "assemblies" / "sample-A"
    assemblies.mkdir(parents=True)
    (assemblies / "assembly.fasta").write_text(
        ">contig 1 circular\nACGTACGT\n>plasmid-2\nTTAA\n", encoding="utf-8"
    )
    split_dir = tmp_path / "contigs"
    merged_dir = tmp_path / "merged"

    subprocess.run(
        [
            sys.executable,
            str(PROJECT / "scripts/02_split_contigs.py"),
            "--assemblies",
            str(tmp_path / "assemblies"),
            "--output",
            str(split_dir),
        ],
        check=True,
    )
    subprocess.run(
        [
            sys.executable,
            str(PROJECT / "scripts/07_merge_contigs.py"),
            "--contigs",
            str(split_dir),
            "--output",
            str(merged_dir),
        ],
        check=True,
    )

    merged = (merged_dir / "sample-A" / "assembly.fasta").read_text(encoding="utf-8")
    assert ">contig 1 circular" in merged
    assert ">plasmid-2" in merged
    assert "ACGTACGT" in merged

