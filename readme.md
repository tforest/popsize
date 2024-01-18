
# popsize module for snpArcher

## Overview

The `popsize` module is an extension of the [snpArcher](https://github.com/harvardinformatics/snpArcher.git) pipeline designed for automated demographic inference, specifically focused on population size change estimation. It streamlines computational tasks, allowing users to concentrate on the biological implications of their analyses and the relevance of the parameters specified using reference tools and methods in this field.

### Purpose

The primary goal of the `popsize` module is to facilitate the estimation of population size changes, addressing key assumptions in population genetics. Users can choose from a suite of powerful demographic inference tools tailored for various data scenarios. The modular design allows customization and parameter adjustments based on specific analysis needs.

### Key Features

1. **Modularity:** The `popsize` module integrates demographic inference tools, offering flexibility and adaptability to different datasets.
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

1. **Installation:**
   - Clone the repository of the `popsize` module under `snpArcher/workflow/modules/`:

     ```bash
     cd snpArcher/workflow/modules/
     git clone https://github.com/tforest/popsize.git
     ```

2. **Configuration:**
   - Adjust settings in the `config.yaml` file in the `config` directory to fit your needs.
   - Important Configurations to Set in `config.yaml`:

     ```yaml
     ## Popsize options
     # popsize modules to use for inference; can be : dadi, swp2, msmc2, psmc, smcpp. Specify multiple tools comma-separated: "dadi, swp2, psmc" for example
     popsize_tools: "dadi"
     # specify generation time for the studied species.
     gen_time: 2
     # specify the average per site mutation rate 
     mut_rate: 3e-8
     # specify whether to use folded SFS or not. 
     folded: True
     # select the number of runs dadi is executing for its optimization phase 
     dadi_optimizations: 100
     ```

3. **Integration with snpArcher:**
   - Add the following lines to the default snpArcher workflow (`snparcher/workflow/Snakefile`):

     ```snakemake
     module popsize:
         snakefile:
             "modules/popsize/Snakefile"
         config:
             config

     use rule * from popsize as popsize_*
     ```

4. **Execution:**
   - Run the snpArcher workflow with the integrated `popsize` module. You can execute on your current snpArcher configuration the `popsize` steps only by running something like:

     ```snakemake
     snakemake --cores 16 --use-conda --conda-frontend mamba --profile ./profiles/slurm --snakefile workflow/Snakefile popsize
     ```

## Dependencies

- Snakemake 7.1x


## Contribution

If you encounter issues or have suggestions, please open an issue on the [GitHub repository](https://github.com/tforest/popsize).

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
