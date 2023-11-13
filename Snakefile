import pandas as pd
import os
include: "common.smk"

#configfile: "config/config.yaml"

samples = pd.read_table(config["samples"], sep=",", dtype=str).replace(' ', '_', regex=True)
REFGENOME = samples['refGenome'].unique().tolist()

rule all:
    input:
        expand("results/{refGenome}/popsize/output_smcpp/{prefix}_model.final.json", refGenome=REFGENOME, prefix=config['final_prefix'])

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
        sfs_file = "results/{refGenome}/popsize/SFS_{prefix}.txt"
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
        sfs_file = "results/{refGenome}/popsize/SFS_{prefix}.txt"
    output:
        stairwayplot2_summary = "results/{refGenome}/popsize/output_stairwayplot2/{prefix}/{prefix}.final.summary"
    conda:
        "envs/deminfhelper.yml"
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --stairwayplot2"

# rule index_vcf:
#     """
#     Index vcf so that it can be used by smcpp
#     """
#     input:
#         vcf = "results/{refGenome}/{prefix}_raw.vcf.gz"
#     output:
#         tabix_index = "results/{refGenome}/{prefix}_raw.vcf.gz.tbi"
#     conda:
#         "envs/deminfhelper.yml"
#     shell:
#         "tabix {input.vcf}"

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
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --smcpp"
