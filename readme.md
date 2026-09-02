# popsize module for snpArcher

## Overview

The `popsize` module is an extension of the [snpArcher](https://github.com/harvardinformatics/snpArcher.git) pipeline designed for automated demographic inference, specifically focused on population size change estimation. It streamlines computational tasks, allowing users to concentrate on the biological implications of their analyses and the relevance of the parameters specified using reference tools and methods in this field.

### Purpose

The primary goal of the `popsize` module is to facilitate the estimation of population size changes, addressing key assumptions in population genetics. Users can choose from a suite of powerful demographic inference tools tailored for various data scenarios. The modular design allows customization and parameter adjustments based on specific analysis needs. In the end, a combined plot presenting an overlay of all the methods used is generated, bringing a comparative visualization to this kind of analysis.

### Key Features

1. **Modularity:** The `popsize` module integrates a selection of population size inference approaches, offering flexibility and adaptability to different datasets and situations.
2. **Transparency:** Parameters used by the integrated tools are transparent, enabling users to understand and control computational aspects.
3. **Statistical Insights:** The module includes essential statistics such as Site Frequency Spectrum (SFS) and Principal Component Analysis (PCA), aiding in population structure evaluation.

### Integrated Tools

1. **PSMC (Pairwise Sequentially Markovian Coalescent)**
   - **Description:** PSMC infers population size history from diploid sequences using a pairwise sequentially Markovian coalescent model. It explores scaled mutation and recombination rates, providing insights into demographic changes.
   - **Reference:** [GitHub - lh3/psmc](https://github.com/lh3/psmc)
   - **Citation:** Li, H., and R. Durbin. 2011. "Inference of human population history from individual whole-genome sequences." Nature 475: 493–496. 

2. **MSMC2 (Multiple Sequentially Markovian Coalescent 2)**
   - **Description:** MSMC2 extends the MSMC model for inferring population size history and separation from whole-genome sequencing data. It provides accurate estimations for a large number of haplotypes.
   - **Reference:** [GitHub - stschiff/msmc2](https://github.com/stschiff/msmc2)
   - **Citation:** Schiffels, S., and K. Wang. 2020. "MSMC and MSMC2: The Multiple Sequentially Markovian Coalescent." In Statistical Population Genomics, edited by J. Y. Dutheil, 147–166. Methods in Molecular Biology, Springer US, New York, NY.

3. **Stairway Plot v2**
   - **Description:** Stairway Plot v2 infers detailed population demographic history using the site frequency spectrum (SFS) from DNA sequence data. It can use both folded and unfolded SFSs and controls for overfitting.
   - **Reference:** [GitHub - xiaoming-liu/stairway-plot-v2](https://github.com/xiaoming-liu/stairway-plot-v2)
   - **Citation:** Liu, X., and Y.-X. Fu, 2020 Stairway Plot 2: demographic history inference with folded SNP frequency spectra. Genome Biology 21: 280.

4. **dadi-cli (Flexible Python Package for Demographic Inference)**
   - **Description:** dadi-cli provides a user-friendly command-line interface for dadi, a flexible Python package for inferring demographic history and the distribution of fitness effects from population genomic data.
   - **Reference:** [dadi-cli Documentation](https://dadi-cli.readthedocs.io)
   - **Citation:** Gutenkunst, R. N., R. D. Hernandez, S. H. Williamson, and C. D. Bustamante. 2009. "Inferring the Joint Demographic History of Multiple Populations from Multidimensional SNP Frequency Data." PLOS Genetics 5: e1000695.
 
5. **SMC++ (Sequentially Markovian Coalescent++)**
   - **Description:** SMC++ is a tool for estimating the size history of populations from whole-genome sequence data. It offers various subcommands for data conversion, estimation, cross-validation, and joint demography analysis.
   - **Reference:** [GitHub - popgenmethods/smcpp](https://github.com/popgenmethods/smcpp)
   - **Citation:** Terhorst, J., J. A. Kamm, and Y. S. Song. 2017. "Robust and scalable inference of population history from hundreds of unphased whole genomes." Nat Genet 49: 303–309.

## Repository Structure

- **bin:** Contains executable files for the demographic inference tools.
- **common.smk:** Shared Snakemake rules used across the pipeline.
- **config:** Configuration files for the module, including `config.yaml` for general settings.
- **envs:** Conda environment files (`deminfhelper.yml`, `smcpp.yml`) specifying dependencies for each step.
- **scripts:** Python scripts (`deminfhelper.py`, `inferences.py`, etc.) for demographic inference and result visualization, including `generate_multihetsep.py` from the `msmc-tools` repository.
- **Snakefile:** The main Snakemake workflow file orchestrating the module.

## Usage

1. **Installation of snpArcher:**
   - This module has been tested with snpArcher v2.2. For now this is the only stable version for which this module is supported.
     ```bash
     wget https://github.com/harvardinformatics/snparcher/archive/refs/tags/v2.2.zip
     # unzip it
     unzip v2.2.zip
     # Rename the directory to make it more convenient
     mv snparcher-2.2 snpArcher
     ```

2. **Configuration of snpArcher:**

- Edit the config/config.yaml file on this lines:

```yaml
samples: "config/samples.csv"
name: "picus_viridis"
source: "GCA_033816785.1"  # Can be a refseq/genbank accession, url, or path
```
 - Create a ```config/samples.csv``` file containing information about your samples, for example:
```csv
sample_id,library_id,input_type,input
SAMN38508702,JF5345,srr,SRR27195338
SAMN38508701,JF5325,srr,SRR27195328
SAMN38508699,JF5258,srr,SRR27195330
SAMN38508698,JF5191,srr,SRR27195331
SAMN38508697,JF5180,srr,SRR27195332
SAMN38508696,MO1993190,srr,SRR27195333
SAMN38508695,MO1991186,srr,SRR27195334
SAMN38508694,MO1991185,srr,SRR27195335
SAMN38508693,MO19711091,srr,SRR27195339
SAMN38508692,MO19711090,srr,SRR27195341
```

*This is an example dataset using Picus viridis samples. (Forest et al., 2024)*

You can find more information about the sample sheet in [snpArcher documentation](https://snparcher.readthedocs.io/en/latest/setup.html#creating-a-sample-sheet).
 - Edit the workflow-profiles/default/config.yaml file on this lines:

```yaml
# slurm
executor: slurm
jobs: 100 # Have up to N jobs submitted at any given time
latency-wait: 900 # Wait N seconds for output files due to latency
retries: 3 # Retry jobs N times.

tmpdir: "/your/tmp/dir"
...
slurm_partition: "YOUR_PARTITION"
...
runtime: 720m # In minutes, here 12h
...
set-resources:
  # Alignment
  bwa_mem: 
    mem_mb: attempt * 16000
    # change runtime to 480m
    runtime: 480m
    
  joint_genomics_db_import:
    mem_mb: attempt * 64000
    mem_mb_reduced: attempt * 48000
    # change runtime to 
    runtime: 1440m
```
     
3. **Do a first run of snpArcher:**

- Create a conda env with snakemake>=9.19.

 - Once you activated the conda environment, while in snpArcher/ top directory, run snpArcher until the end of the process. (ie. creation of the VCFs):

```bash
snakemake --profile ./workflow-profiles/default
```

2. **Installation of the *popsize* module:**	
- Clone the repository of the `popsize` module under `snpArcher/workflow/modules/`:

```bash
cd snpArcher/workflow/modules/
git clone https://github.com/tforest/popsize.git
```

3. **Integration with snpArcher:**

- In the snpArcher general Snakefile (`snparcher/workflow/Snakefile`), edit  ```rule all``` under the "Default target" category :
```yaml
rule all:
    """Default target: full pipeline."""
    default_target: True
    input:
        RAW_VCF,
        ...
        # ADD THIS LINE
        "results/popsize/.done",
```

- And add the following lines to the default snpArcher workflow (`snparcher/workflow/Snakefile`):

```yaml
module popsize:
snakefile:
        "modules/popsize/Snakefile"
    config:
        config

use rule * from popsize exclude all as popsize_*
```

4. **Configuration of the popsize module:**
   - Adjust settings in the `snpArcher/workflow/modules/popsize/config/config.yaml` file to fit your needs.
   - Here is an example `config.yaml`:

```yaml
### POP. DEFINITION
pop_name: picus_viridis

# generation time (avg. reprod. time in years)
gen_time: 5.6
# mutation rate (per base per generation)
mut_rate: 5e-09

### FILTERS
# contigs longer than this are kept
length_cutoff: 10000
# regex; keep only contigs whose name matches
contig_filter: "CM.*"

### DEFAULT RESOURCES
mem: 4096

### PLOTTING
# applies to every plot except PSMC's, which is rendered by psmc's own
# perl script (plot_psmc.pl) and is always .eps
plot_format: svg

### SFS
folded: True

### PCA
n_clust_kmeans: 3

### TOOLS
# Inferences to run (full list or subset): dadi, swp2, psmc, msmc2, smcpp
popsize_tools: dadi,swp2,psmc,msmc2,smcpp

# Dadi
optimizations: 1000
p0: 1, 1, 0.2, 1
lower_bound: 0.01, 0.01, 0.005, 0.1
upper_bound: 10, 4, 5, 10

# MSMC2
msmc2_kwargs: -i 25 -p 1*2+25*1+1*2+1*3

# PSMC
psmc_kwargs: -N25 -t15 -r5 -p "4+25*2+4+6"
# -X caps the x-axis ()
# -Y caps the y-axis (effective pop size, x1e4).
plot_psmc_kwargs: -x 10**4 -Y 11

# Stairway Plot 2
path_to_stairwayplot2: workflow/modules/popsize/bin/stairway_plot_es/
blueprint_template: workflow/modules/popsize/config/swp2_template.blueprint
```

- Once snpArcher executed properly and all the files are created (specifically the VCFs under results/vcfs), you can run popsize by re-running snpArcher and it will generate the missing output created by the module:

```bash
snakemake --profile ./workflow-profiles/default
```

## Dependencies

- Snakemake >= 9.19


## Contribution

If you encounter issues or have suggestions, please open an issue on the [GitHub repository](https://github.com/tforest/popsize/issues).

## Acknowledgments

The `popsize` module integrates scripts from the [deminfhelper project](https://github.com/tforest/deminfhelper.git) and [msmc-tools project](https://github.com/stschiff/msmc-tools.git). Please visit [snpArcher Github repository](https://github.com/harvardinformatics/snpArcher.git) for more details on the robust framework that forms the basis of this module.

## Citation

If you use the snpArcher, please consider citing the following publication:

- Mirchandani, C. D., A. J. Shultz, G. W. C. Thomas, S. J. Smith, M. Baylis et al., 2023. "A fast, reproducible, high-throughput variant calling workflow for population genomics." Molecular Biology and Evolution, msad270.

Popsize integrated tools rely on these published works:

- Gutenkunst, R. N., R. D. Hernandez, S. H. Williamson, and C. D. Bustamante. 2009. "Inferring the Joint Demographic History of Multiple Populations from Multidimensional SNP Frequency Data." PLOS Genetics 5: e1000695.
- Li, H., and R. Durbin. 2011. "Inference of human population history from individual whole-genome sequences." Nature 475: 493–496.
- Schiffels, S., and K. Wang. 2020. "MSMC and MSMC2: The Multiple Sequentially Markovian Coalescent." In Statistical Population Genomics, edited by J. Y. Dutheil, 147–166. Methods in Molecular Biology, Springer US, New York, NY.
- Terhorst, J., J. A. Kamm, and Y. S. Song. 2017. "Robust and scalable inference of population history from hundreds of unphased whole genomes." Nat Genet 49: 303–309.
- Liu, X., and Y.-X. Fu, 2020 Stairway Plot 2: demographic history inference with folded SNP frequency spectra. Genome Biology 21: 280.

The example dataset is taken from this study:
- Forest, T., Achaz, G., Marbouty, M., Bignaud, A., Thierry, A., Koszul, R., Milhes, M., Lledo, J., Pons, J.M., & Fuchs, J. (2024). Chromosome-level genome assembly of the European green woodpecker Picus viridis. _G3: Genes, Genomes, Genetics_, _14_(5), jkae042.
