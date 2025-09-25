import sys


def test_cli_main_uses_provided_config(monkeypatch, tmp_path):
    # Create a dummy config path (not used by fake pipeline)
    cfg = tmp_path / "config.toml"
    cfg.write_text("log_level='INFO'\nlog_file='x.log'\n")

    # Prepare capture variables
    called = {}

    class FakePipeline:
        def __init__(self, config_path: str):
            called["config_path"] = config_path

        def run(self):
            called["ran"] = True

    import foldfusion.foldfusion as ff

    monkeypatch.setattr(ff, "FoldFusion", FakePipeline)
    monkeypatch.setattr(sys, "argv", ["foldfusion", str(cfg)])

    # Invoke
    ff.main()

    assert called.get("ran") is True
    assert called.get("config_path") == str(cfg)
