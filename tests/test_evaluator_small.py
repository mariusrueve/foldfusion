import json
from pathlib import Path

from foldfusion.evaluation.evaluator import Evaluator


def _pdb_atom_line(n, name, res, chain, resno, x, y, z, elem):
    # Construct a minimal fixed-width PDB ATOM record matching parser slices
    return (
        f"ATOM  {n:5d} {name:<4}{res:>3} {chain:1}{resno:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{50.00:6.2f}          {elem:>2}\n"
    )


def _write_minimal_pdb(path: Path, shift: float = 0.0):
    lines = [
        _pdb_atom_line(1, "CA", "ALA", "A", 1, 0.000 + shift, 0.000, 0.000, "C"),
        _pdb_atom_line(2, "CA", "GLY", "A", 2, 1.000 + shift, 0.000, 0.000, "C"),
    ]
    path.write_text("".join(lines))


def _write_minimal_sdf(path: Path):
    # Simple 1-atom molecule at origin
    header = [
        "TestMol\n",
        "Created by test\n",
        "\n",
    ]
    counts = f"{1:>3}{0:>3}  0  0  0  0            999 V2000\n"
    # atom line: x y z and element at columns 31:34
    atom = f"{0.0:>10.4f}{0.0:>10.4f}{0.0:>10.4f}  C  0  0  0  0  0  0  0  0  0  0\n"
    path.write_text("".join(header + [counts, atom, "M  END\n$$$$\n"]))


def test_evaluator_writes_metrics(tmp_path):
    out_dir = tmp_path / "out"
    af = tmp_path / "af.pdb"
    exp = tmp_path / "exp.pdb"
    lig = tmp_path / "lig.sdf"
    _write_minimal_pdb(af, shift=0.0)
    _write_minimal_pdb(exp, shift=0.1)
    _write_minimal_sdf(lig)

    ev = Evaluator(out_dir)
    siena_structs = [
        {
            "pdb_code": "XXXX",
            "chain": "A",
            "ensemble_path": exp,
            "ligand_pdb_code": "LIG",
            "backbone_rmsd": 0.1,
            "all_atom_rmsd": 0.2,
        }
    ]
    ligand_structs = {"XXXX": [{"ligand_id": "LIG_A_1", "path": str(lig)}]}

    ev.evaluate("QTEST", "pre-jamda", af, siena_structs, ligand_structs)

    eval_json = ev.output_dir / "evaluation.json"
    assert eval_json.exists()
    data = json.loads(eval_json.read_text())
    assert "QTEST" in data
    assert "XXXX" in data["QTEST"]
    # Confirm metric keys present
    lig_map = data["QTEST"]["XXXX"]["LIG_A_1"]["pre-jamda"]
    assert "local_rmsd" in lig_map
    assert "tcs" in lig_map
