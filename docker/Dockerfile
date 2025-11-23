# Start from the official Ubuntu 24.04.2 LTS base image
FROM ubuntu:24.04

# Set environment variables to avoid interactive prompts during install
ENV DEBIAN_FRONTEND=noninteractive

# Install essential packages
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    vim \
    sudo \
    pip \
    build-essential \
    python3-pandas \
    pkg-config \
    zip \
    cmake \
    libprotobuf-dev \
    protobuf-compiler \
    libtbb-dev \
    python3-pip \
    git \
    nodejs \
    pigz \
    minimap2 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Download and install Miniforge (Conda)
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" \
    && bash Miniforge3-$(uname)-$(uname -m).sh -b \
    && rm Miniforge3-$(uname)-$(uname -m).sh

# Set Conda/Miniforge paths
ENV PATH=/root/miniforge3/bin:$PATH

# Initialize Conda and ensure it’s sourced on login
RUN /root/miniforge3/bin/conda init bash && \
    echo "source ~/.bashrc" >> /root/.bash_profile

# Install Snakemake, taxoniumtools, and matplotlib using pip
RUN pip install snakemake taxoniumtools matplotlib

# Install viral_usher using Conda
RUN conda install -c bioconda -y viral-usher

# Create workspace directory
RUN mkdir -p /workspace

# Clone metaWEPP repository
RUN git clone https://github.com/TurakhiaLab/metaWEPP.git

# Clone WEPP repository with submodules
WORKDIR /metaWEPP
RUN git clone --recurse-submodules https://github.com/TurakhiaLab/WEPP.git

# Clone and install Kraken2
RUN git clone https://github.com/DerrickWood/kraken2.git && \
    cd kraken2 && \
    ./install_kraken2.sh . && \
    echo 'export PATH="/metaWEPP/kraken2:$PATH"' >> /root/.bashrc


# Default command: start an interactive login shell to source ~/.bashrc
CMD [ "bash", "--login" ]
