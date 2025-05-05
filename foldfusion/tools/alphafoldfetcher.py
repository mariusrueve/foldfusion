from foldfusion.tools.tool import Tool
import requests
from pathlib import Path


class AlphaFoldFetcher(Tool):
    def __init__(self, config):
        super().__init__(config)
        self.tool_name = "alphafoldfetcher"
        self.uniprot_id = self.config.get("run_settings", "")["uniprot_id"]
        self.output_dir = self.output_dir / "alphafold"

    def get_alphafold_model(self) -> Path:
        self.output_dir.mkdir(parents=True, exist_ok=True)

        versions = ["v4", "v3"]
        last_error = None

        for version in versions:
            model_version_suffix = f"model_{version}"
            url = (
                f"https://alphafold.ebi.ac.uk/files/"
                f"AF-{self.uniprot_id}-F1-{model_version_suffix}.pdb"
            )
            output_file_name = f"AF-{self.uniprot_id}-F1-{model_version_suffix}.pdb"
            output_file_path = self.output_dir / output_file_name

            try:
                response = requests.get(url, timeout=60)
                response.raise_for_status()

                if response.status_code == 200:
                    output_file_path.write_text(response.text, encoding="utf-8")
                    pdb_path = output_file_path
                    last_error = None
                    break
            except requests.exceptions.HTTPError as e:
                last_error = e
                if e.response.status_code == 404:
                    print(
                        f"Model version {version} not found for ",
                        f"{self.uniprot_id} (404)",
                    )
                else:
                    print(
                        f"HTTP error fetching structure version {version} for "
                        f"{self.uniprot_id}. Status Code: {e.response.status_code}. "
                        f"Error: {e}"
                    )
            except requests.exceptions.RequestException as e:
                last_error = e
                print(
                    f"Failed to fetch structure version {version} for "
                    f"{self.uniprot_id} due to a network or request issue. "
                    f"Error: {e}"
                )
        if pdb_path:
            return pdb_path.resolve()
        else:
            error_message = (
                f"Failed to fetch any structure version for {self.uniprot_id} "
                f"after trying {versions}."
            )
            print(error_message)
            if last_error:
                print(f"Last error encountered: {last_error}")
                raise FileNotFoundError(
                    f"{error_message} Last error: {last_error}"
                ) from last_error
            else:
                raise FileNotFoundError(error_message)
