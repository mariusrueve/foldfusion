from foldfusion.tools.dogsite3 import Dogsite3
from foldfusion.tools.siena import Siena
from foldfusion.tools.siena_db import SienaDB


def test_dogsite3_command(tmp_path):
    exe = tmp_path / "dogsite3"
    exe.write_text("")
    pdb = tmp_path / "model.pdb"
    pdb.write_text(
        "ATOM      1  C   ALA A   1      0.000   0.000   0.000  1.00 50.00"
        "           C\n"
    )
    out = tmp_path / "out"
    d = Dogsite3(exe, pdb, out)
    cmd = d.get_command()
    # Command begins with executable, contains proteinFile arg and writeSiteResiduesEDF
    assert str(exe) in cmd[0]
    assert "--proteinFile" in cmd
    assert "--writeSiteResiduesEDF" in cmd


def test_siena_command(tmp_path):
    exe = tmp_path / "siena"
    exe.write_text("")
    edf = tmp_path / "out.edf"
    edf.write_text("REFERENCE X\n")
    db = tmp_path / "siena_db"
    db.write_text("x" * 2048)  # >1KB considered valid
    pdb_dir = tmp_path / "pdb"
    pdb_dir.mkdir()
    out = tmp_path / "out"
    s = Siena(exe, edf, db, pdb_dir, out)
    cmd = s.get_command()
    # Check required flags present
    assert cmd[0].endswith("siena")
    assert "--edf" in cmd and str(edf.resolve()) in cmd
    assert "--database" in cmd and str(db.resolve()) in cmd
    assert "--output" in cmd


def test_siena_db_command_and_naming(tmp_path):
    exe = tmp_path / "gendb"
    exe.write_text("")
    pdb_dir = tmp_path / "pdb"
    pdb_dir.mkdir()
    out = tmp_path / "out"
    # Case 1: configured with .db suffix
    db_path = out / "mydb.db"
    sdb = SienaDB(exe, db_path, pdb_dir, 1, out)
    cmd = sdb.get_command()
    # Ensure the database argument uses the provided name including suffix
    assert "--database" in cmd
    assert cmd[cmd.index("--database") + 1] == db_path.name
    # Case 2: no suffix
    db_path2 = out / "mydb"
    sdb2 = SienaDB(exe, db_path2, pdb_dir, 1, out)
    cmd2 = sdb2.get_command()
    assert cmd2[cmd2.index("--database") + 1] == db_path2.name
