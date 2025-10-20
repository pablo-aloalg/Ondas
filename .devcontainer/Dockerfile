# Base Dev Container with Miniconda
FROM mcr.microsoft.com/devcontainers/miniconda:1-3

# Copy environment.yml into the container
COPY environment.yml /tmp/environment.yml

# Ensure Conda environment auto-activates in every terminal
RUN echo "source /opt/conda/etc/profile.d/conda.sh && conda activate base" >> /home/vscode/.bashrc
