# Ensure config.toml exists and contains the correct paths/options
./main.py \
  --uniprot_id Q9Y233 \
  --siena_db siena_db.db \
  --output_dir pipeline_output/ \
  --config config.toml \
  # Add -v for verbose logging if needed
  # -v