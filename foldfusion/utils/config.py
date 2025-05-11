from pathlib import Path
import tomllib


class Config:
    # TODO: Do more in the class itself and offload here from tool classes
    def __init__(self, config_path: Path):
        self.config_path = config_path
        self.dict = None
        self.load_config()

    def load_config(self):
        if not self.config_path.exists():
            raise FileNotFoundError(f"Config file not found at {self.config_path}")

        try:
            with open(
                self.config_path, "rb"
            ) as f:  # TOML must be opened in binary mode
                config = tomllib.load(f)
            self.dict = config
        except Exception as e:
            raise ValueError(f"Failed to parse config file: {e}")

    def get_tool_config(self, tool_name: str):
        return self.dict.get(tool_name, {})


if __name__ == "__main__":
    c = Config(Path("config.toml"))
    print(c.get_tool_config("siena_db"))
