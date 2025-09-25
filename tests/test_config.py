from pathlib import Path

import pytest

from foldfusion.utils.config import Config


def write_minimal_valid_config(tmp_path: Path) -> Path:
    # Create dummy executables and directories
    exes = {}
    for name in [
        "dogsite3",
        "siena",
        "siena_db",
        "ligand_extractor",
        "jamda_scorer",
    ]:
        p = tmp_path / f"{name}.bin"
        p.write_text("#!/bin/sh\n")
        exes[name] = p

    pdb_dir = tmp_path / "pdb"
    pdb_dir.mkdir()

    cfg = tmp_path / "config.toml"
    cfg.write_text(
        "\n".join(
            [
                'log_level = "INFO"',
                'log_file = "ff.log"',
                'uniprot_ids = ["Q8CA95"]',
                f'output_dir = "{(tmp_path / "out").as_posix()}"',
                "pipeline_concurrency = 1",
                f'dogsite3_executable = "{exes["dogsite3"].as_posix()}"',
                f'siena_executable = "{exes["siena"].as_posix()}"',
                f'siena_db_executable = "{exes["siena_db"].as_posix()}"',
                # Use default DB location; ensure pdb_directory exists
                f'pdb_directory = "{pdb_dir.as_posix()}"',
                "pdb_format = 1",
                (
                    "ligand_extractor_executable = "
                    f'"{exes["ligand_extractor"].as_posix()}"'
                ),
                (f'jamda_scorer_executable = "{exes["jamda_scorer"].as_posix()}"'),
            ]
        )
    )
    return cfg


def test_config_loads_ids_and_paths(tmp_path):
    cfg_path = write_minimal_valid_config(tmp_path)
    c = Config(cfg_path)

    assert c.uniprot_ids == ["Q8CA95"]
    assert c.output_dir.exists()
    assert c.pipeline_concurrency == 1
    # executables should resolve and exist
    assert c.dogsite3_executable.exists()
    assert c.siena_executable.exists()
    assert c.siena_db_executable.exists()
    assert c.ligand_extractor_executable.exists()
    assert c.jamda_scorer_executable.exists()
    # pdb directory and format
    assert c.pdb_directory.exists()
    assert c.pdb_format in (0, 1)


def test_config_rejects_both_ids_sources(tmp_path):
    cfg = tmp_path / "config.toml"
    ids_file = tmp_path / "ids.txt"
    ids_file.write_text("Q9Y233\n")
    cfg.write_text(
        "\n".join(
            [
                'log_level = "INFO"',
                'log_file = "ff.log"',
                'uniprot_ids = ["Q8CA95"]',
                f'uniprot_ids_file = "{ids_file.as_posix()}"',
                f'output_dir = "{(tmp_path / "out").as_posix()}"',
                # Minimal keys to trigger error before executable checks
                f'pdb_directory = "{(tmp_path / "pdb").as_posix()}"',
                "pdb_format = 1",
            ]
        )
    )

    with pytest.raises(ValueError):
        Config(cfg)


def test_config_ids_from_file(tmp_path):
    ids_file = tmp_path / "ids.txt"
    ids_file.write_text("Q9Y233\n\nQ9QYJ6\n")
    cfg = tmp_path / "config.toml"
    pdb_dir = tmp_path / "pdb"
    pdb_dir.mkdir()
    cfg.write_text(
        "\n".join(
            [
                'log_level = "INFO"',
                'log_file = "ff.log"',
                f'uniprot_ids_file = "{ids_file.as_posix()}"',
                f'output_dir = "{(tmp_path / "out").as_posix()}"',
                f'pdb_directory = "{pdb_dir.as_posix()}"',
                "pdb_format = 1",
            ]
        )
    )
    c = Config(cfg)
    assert c.uniprot_ids == ["Q9Y233", "Q9QYJ6"]
