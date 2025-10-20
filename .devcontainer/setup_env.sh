#!/bin/bash

# Activate Conda
source /opt/conda/etc/profile.d/conda.sh

# Create or update environment from environment.yml
conda env create -f /tmp/environment.yml || conda env update -f /tmp/environment.yml --prune

# Auto-activate the environment in all new terminals
ENV_NAME=$(head -1 /tmp/environment.yml | cut -d' ' -f2)
echo "source /opt/conda/etc/profile.d/conda.sh && conda activate $ENV_NAME" >> ~/.bashrc

# Trust all notebooks in the repo
jupyter trust $(find . -name '*.ipynb')

echo "✅ Conda environment installed and all notebooks trusted!"
