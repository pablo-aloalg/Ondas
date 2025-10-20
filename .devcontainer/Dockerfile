# Use the official Miniconda Dev Container image from Microsoft
FROM mcr.microsoft.com/devcontainers/miniconda:1-3

# Copy environment.yml into a temp folder
COPY environment.yml /tmp/conda-tmp/environment.yml

# Create/Update the 'base' environment (can change 'base' to another name if desired)
RUN umask 0002 && \
    source /opt/conda/etc/profile.d/conda.sh && \
    conda env update -n base -f /tmp/conda-tmp/environment.yml && \
    conda clean -afy && \
    rm -rf /tmp/conda-tmp

# Automatically activate Conda in all new terminals
RUN echo "source /opt/conda/etc/profile.d/conda.sh && conda activate base" >> /home/vscode/.bashrc
