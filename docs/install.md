# Install

metaWEPP offers multiple installation methods.

1. Bioconda (Recommended)
2. Dockerfile
3. Shell Commands

## **Option-1: Install via Bioconda (Recommended)** <a name=conda></a>

**Step 1:** Create a new conda environment for metaWEPP.
```bash
conda create --name metawepp-env
conda activate metawepp-env
conda config --env --add channels bioconda
conda config --env --add channels conda-forge
conda config --env --set channel_priority flexible
conda install metawepp
```

!!!Note
    ⚠️ You can use `conda install metawepp --solver=libmamba` to enable a faster dependency resolution and installation.

**Step 2:** Confirm proper working by running the following command. This should print metaWEPP's help menu.
```bash
run-metawepp help --cores 1
```

**Step 3:** Create a `data` directory and start analyzing your samples with metaWEPP. If you are running samples from multiple `data` directories, specify the `.snakemake` directory created in one run as the `--conda-prefix` for the others to avoid redundant creation of Snakemake conda environments.

Before trying the [examples](quickstart.md#example), please ensure that you have downloaded the `simulated_metagenomic_sample` into the `data` directory from [here](https://github.com/TurakhiaLab/metaWEPP/tree/main/data/simulated_metagenomic_sample).

!!!Note
    ⚠️ If you plan to generate a MAT for any species within the metaWEPP workflow using viral_usher, Docker access must be available on your system.


## **Option-2: Install via Dockerfile** <a name=docker></a>

**Step 1:** Clone the metaWEPP repository.
```bash
git clone https://github.com/TurakhiaLab/metaWEPP.git
cd metaWEPP
```

**Step 2:** Build a Docker Image.
```bash
cd docker
docker build -t metawepp .
cd ..
```

**Step 3:** Start and run Docker container. The command below will take you inside the Docker container with metaWEPP already installed.
```bash
# -p <host_port>:<container_port> → Maps container port to a port on your host (Accessing Dashboard, NOT needed otherwise)
# Replace <host_port> with your desired local port (e.g., 100 or 8080)
docker run -it -p 80:80 metawepp
```

**Step 4:** Confirm proper working by running the following command. This should print metaWEPP's help menu.
```bash
run-metawepp help --cores 1
```

All set to try the [examples](quickstart.md#example).

!!!Note
    ⚠️ If you do not already have a MAT for a pathogen and wish to generate one using viral_usher, you must run viral_usher outside the metaWEPP Docker container and then copy the resulting MAT into the container. This is necessary because viral_usher invokes its own Docker instance during execution and will fail when run from within another Docker container.


## **Option-3: Install via Shell Commands (requires sudo access)** <a name=script></a>

**Step 1:** Clone the repository.
```bash
git clone https://github.com/TurakhiaLab/metaWEPP.git
cd metaWEPP
chmod +x run-metawepp
```

**Step 2:** Update `~/.bashrc` for linux or `~/.zshrc` for macOS.
```bash
echo "
run-metawepp() {
    snakemake -s $PWD/Snakefile \"\$@\"
}
export -f run-metawepp
" >> ~/.bashrc

source ~/.bashrc
```

**Step 3:** Install Kraken.
The following commands install kraken and also update the `$PATH` variable for running the tool easily.
```bash
git clone https://github.com/DerrickWood/kraken2.git
cd kraken2
./install_kraken2.sh .
echo -e "\nexport PATH=\"$(pwd):\$PATH\"" >> ~/.bashrc
source ~/.bashrc
cd ..
```

**Step 4:** Install `Minimap2`, `viral_usher`, `matplotlib`, and `snakemake`.
```bash
sudo apt-get install minimap2
pip install viral_usher matplotlib snakemake
```

**Step 5:** Install `WEPP`.
```bash
git clone --recurse-submodules https://github.com/TurakhiaLab/WEPP.git
cd WEPP
chmod +x run-wepp
```

View the WEPP installation guide starting from Option 3 of the [WEPP](https://github.com/TurakhiaLab/WEPP/tree/main?tab=readme-ov-file#-option-3-install-via-shell-commands-requires-sudo-access) repository.

**Step 6:** Confirm proper working by running the following command. This should print metaWEPP's help menu.
```bash
run-metawepp help --cores 1
```

All set to try the [examples](quickstart.md#example).

!!!Note
    ⚠️ If you plan to generate a MAT for any species within the metaWEPP workflow using viral_usher, Docker access must be available on your system.
