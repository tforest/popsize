import pandas as pd
import os
import yaml

# Include additionnal functions
include: "common.smk"

# Add complementary config
popsize_config = yaml.safe_load(open("workflow/modules/popsize/config/config.yaml"))
config.update(popsize_config)

samples = pd.read_table(config["samples"], sep=",", dtype=str).replace(' ', '_', regex=True)
REFGENOME = samples['refGenome'].unique().tolist()


rule all:
    input:
        expand("results/{refGenome}/popsize/SFS_{prefix}.fs", refGenome=REFGENOME, prefix=config['final_prefix']),
        expand("results/{refGenome}/popsize/output_stats/{prefix}_PCA.png", refGenome=REFGENOME, prefix=config['final_prefix']),
        *(expand("results/{refGenome}/popsize/output_smcpp/{prefix}_model.final.json", refGenome=REFGENOME, 
        prefix=config['final_prefix']) if "smcpp" in config['popsize_tools'] else []),
        *(expand("results/{refGenome}/popsize/output_stairwayplot2/{prefix}/{prefix}.final.summary", refGenome=REFGENOME, 
        prefix=config['final_prefix']) if "swp2" in config['popsize_tools'] else []),
        # Wait for dadi output if used in config
        *(expand("results/{refGenome}/popsize/output_dadi/{prefix}.InferDM.bestfits", refGenome=REFGENOME, 
        prefix=config['final_prefix']) if "dadi" in config['popsize_tools'] else []),
        # Wait for PSMC output if used in config
        *(expand("results/{refGenome}/popsize/output_psmc/{prefix}_combined.psmc.final", refGenome=REFGENOME, 
        prefix=config['final_prefix']) if "psmc" in config['popsize_tools'] else []),
        # Wait for MSMC2 output if used in config
        *(expand("results/{refGenome}/popsize/output_msmc2/{prefix}_msmc2.final.txt", refGenome=REFGENOME, 
        prefix=config['final_prefix']) if "msmc2" in config['popsize_tools'] else [])

rule init_module:
    """
    intial checks and creates the config file
    """
    input:
        vcf = "results/{refGenome}/{prefix}_raw.vcf.gz",
        results_folder = "results/{refGenome}",
        config_template = "workflow/modules/popsize/config/deminfhelper_template.yml"
    output:
        config_file = "results/{refGenome}/popsize/{prefix}_deminfhelper.yml"
    run:
        build_deminfhelper_config(results_folder = input.results_folder, vcf_file = input.vcf,
        prefix = config['final_prefix'], output = output.config_file, global_config = config, module_config = input.config_template)

rule compute_sfs:
    """
    computes sfs from the VCF
    """
    input:
        config_file = "results/{refGenome}/popsize/{prefix}_deminfhelper.yml"
    output:
        sfs_file = "results/{refGenome}/popsize/SFS_{prefix}.fs"
    conda:
        "envs/deminfhelper.yml"
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --sfs --config_file {input.config_file}"

rule stairwayplot2:
    """
    run stairwayplot2 using the generated config
    """
    input:
        config_file = "results/{refGenome}/popsize/{prefix}_deminfhelper.yml",
        sfs_file = "results/{refGenome}/popsize/SFS_{prefix}.fs"
    output:
        stairwayplot2_summary = "results/{refGenome}/popsize/output_stairwayplot2/{prefix}/{prefix}.final.summary"
    conda:
        "envs/deminfhelper.yml"
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --stairwayplot2 ;"+ \
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --plot_stairwayplot2" 

rule smcpp:
    """
    run smc++ using the generated config
    """
    input:
        config_file = "results/{refGenome}/popsize/{prefix}_deminfhelper.yml",
        vcf_index = "results/{refGenome}/{prefix}_raw.vcf.gz.tbi"
    output:
        smcpp_summary = "results/{refGenome}/popsize/output_smcpp/{prefix}_model.final.json"
    conda:
        "envs/smcpp.yml"
    # TODO: CHECK how to make the RAM scale up
    resources: mem_mb=10000, mem_mib=1908

    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --smcpp;"+ \
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --plot_smcpp"
        
rule dadi:
    """
    run dadi using the generated config
    """
    input:
        config_file = "results/{refGenome}/popsize/{prefix}_deminfhelper.yml",
        sfs_file = "results/{refGenome}/popsize/SFS_{prefix}.fs"
    output:
        dadi_summary = "results/{refGenome}/popsize/output_dadi/{prefix}.InferDM.bestfits"
    conda:
        "envs/deminfhelper.yml"
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --dadi ;"+ \
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --plot_dadi"

rule psmc:
    """
    run psmc using the generated config
    """
    input:
        config_file = "results/{refGenome}/popsize/{prefix}_deminfhelper.yml",
        vcf = "results/{refGenome}/{prefix}_raw.vcf.gz",
        vcf_index = "results/{refGenome}/{prefix}_raw.vcf.gz.tbi",
        ref_genome = "results/{refGenome}/data/genome/{refGenome}.fna"
    output:
        psmc_output = "results/{refGenome}/popsize/output_psmc/{prefix}_combined.psmc.final"
    conda:
        "envs/deminfhelper.yml"
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --psmc ;"+ \
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --plot_psmc"

rule msmc2:
    """
    run msmc2 using the generated config
    """
    input:
        config_file = "results/{refGenome}/popsize/{prefix}_deminfhelper.yml",
        vcf = "results/{refGenome}/{prefix}_raw.vcf.gz"
    output:
        msmc2_output = "results/{refGenome}/popsize/output_msmc2/{prefix}_msmc2.final.txt"
    conda:
        "envs/deminfhelper.yml"
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --msmc2 ;"+ \
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --plot_msmc2"

rule statistics:
    """
    Compute some statistics like PCA or Genotyping Quality distribution
    """
    input:
        config_file = "results/{refGenome}/popsize/{prefix}_deminfhelper.yml",
        vcf = "results/{refGenome}/{prefix}_raw.vcf.gz"
    output:
        pca_output = "results/{refGenome}/popsize/output_stats/{prefix}_PCA.png"
    conda:
        "envs/deminfhelper.yml"
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --pca ;"+ \
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --gq_distrib"


        
