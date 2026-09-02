include: "common.smk"
include: "config/resources.smk"

PARAMS = resolve_popsize_params()
PREFIX, _TOOLS = get_popsize_targets(PARAMS)
VCF, REF_FASTA = get_snparcher_paths(config)

_ALL_TARGETS = [
    f"results/popsize/SFS_{PREFIX}.fs",
    f"results/popsize/output_stats/{PREFIX}_PCA.{PARAMS['plot_format']}",
    f"results/popsize/output_smcpp/{PREFIX}.final.json" if "smcpp" in _TOOLS else [],
    f"results/popsize/output_stairwayplot2/{PREFIX}/{PREFIX}.final.summary" if "swp2" in _TOOLS else [],
    f"results/popsize/output_dadi/{PREFIX}.InferDM.bestfits" if "dadi" in _TOOLS else [],
    f"results/popsize/output_psmc/{PREFIX}.eps" if "psmc" in _TOOLS else [],
    f"results/popsize/output_msmc2/{PREFIX}_msmc2.final.txt" if "msmc2" in _TOOLS else [],
    f"results/popsize/{PREFIX}_combined_plot.{PARAMS['plot_format']}",
]

rule all:
    input:
        _ALL_TARGETS


rule done:
    input:
        _ALL_TARGETS
    output:
        touch("results/popsize/.done"),


rule decompress_reference:
    input:
        ref=REF_FASTA,
    output:
        ref="results/popsize/ref.fa",
    shell:
        "gzip -dc {input.ref} > {output.ref}"


rule init_module:
    input:
        vcf=VCF,
    output:
        config_file=f"results/popsize/{PREFIX}_deminfhelper.yml",
    run:
        build_deminfhelper_config(
            out_dir="results/popsize",
            vcf_file=input.vcf,
            ref_fasta="results/popsize/ref.fa",
            prefix=PREFIX,
            output=output.config_file,
            params=PARAMS,
        )


rule compute_sfs:
    input:
        config_file=ancient(f"results/popsize/{PREFIX}_deminfhelper.yml"),
    output:
        sfs_file=f"results/popsize/SFS_{PREFIX}.fs",
    conda:
        "envs/deminfhelper.yml"
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --sfs --config_file {input.config_file} ;"
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --plot_sfs"


rule stairwayplot2:
    input:
        config_file=ancient(f"results/popsize/{PREFIX}_deminfhelper.yml"),
        sfs_file=f"results/popsize/SFS_{PREFIX}.fs",
    output:
        stairwayplot2_summary=f"results/popsize/output_stairwayplot2/{PREFIX}/{PREFIX}.final.summary",
    conda:
        "envs/deminfhelper.yml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * resources['swp2']['mem_mb']
    threads:
        resources['swp2']['threads']
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --stairwayplot2 ;"
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --plot_stairwayplot2"


rule smcpp:
    input:
        config_file=ancient(f"results/popsize/{PREFIX}_deminfhelper.yml"),
        vcf_index=f"{VCF}.tbi",
    output:
        smcpp_summary=f"results/popsize/output_smcpp/{PREFIX}.final.json",
    conda:
        "envs/smcpp.yml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * resources['smcpp']['mem_mb']
    threads:
        resources['smcpp']['threads']
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --cpus {threads} --smcpp ;"
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --plot_smcpp"


rule dadi:
    input:
        config_file=ancient(f"results/popsize/{PREFIX}_deminfhelper.yml"),
        sfs_file=f"results/popsize/SFS_{PREFIX}.fs",
    output:
        dadi_summary=f"results/popsize/output_dadi/{PREFIX}.InferDM.bestfits",
    conda:
        "envs/deminfhelper.yml"
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --dadi ;"
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --plot_dadi"


rule psmc:
    input:
        config_file=ancient(f"results/popsize/{PREFIX}_deminfhelper.yml"),
        vcf=VCF,
        vcf_index=f"{VCF}.tbi",
        ref_genome="results/popsize/ref.fa",
    output:
        psmc_output=f"results/popsize/output_psmc/{PREFIX}_combined.psmc.final",
    params:
        psmc_kwargs=PARAMS["psmc_kwargs"],
    resources:
        mem_mb=lambda wildcards, attempt: attempt * resources['psmc']['mem_mb']
    threads:
        resources['psmc']['threads']
    conda:
        "envs/deminfhelper.yml"
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} "
        "--psmc_kwargs '{params.psmc_kwargs}' --psmc"


rule psmc_plot:
    input:
        config_file=ancient(f"results/popsize/{PREFIX}_deminfhelper.yml"),
        psmc_output=f"results/popsize/output_psmc/{PREFIX}_combined.psmc.final",
    output:
        psmc_figure=f"results/popsize/output_psmc/{PREFIX}.eps",
    params:
        plot_psmc_kwargs=PARAMS["plot_psmc_kwargs"],
    conda:
        "envs/deminfhelper.yml"
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} "
        "--plot_psmc_kwargs '{params.plot_psmc_kwargs}' --plot_psmc"


rule combined_plot:
    input:
        config_file=ancient(f"results/popsize/{PREFIX}_deminfhelper.yml"),
        dadi=f"results/popsize/output_dadi/{PREFIX}.InferDM.bestfits" if "dadi" in _TOOLS else [],
        smcpp=f"results/popsize/output_smcpp/{PREFIX}.final.json" if "smcpp" in _TOOLS else [],
        swp2=f"results/popsize/output_stairwayplot2/{PREFIX}/{PREFIX}.final.summary" if "swp2" in _TOOLS else [],
        msmc2=f"results/popsize/output_msmc2/{PREFIX}_msmc2.final.txt" if "msmc2" in _TOOLS else [],
        psmc=f"results/popsize/output_psmc/{PREFIX}.eps" if "psmc" in _TOOLS else [],
    output:
        combined_plot=f"results/popsize/{PREFIX}_combined_plot.{PARAMS['plot_format']}",
    conda:
        "envs/deminfhelper.yml"
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --combined_plot"


rule msmc2:
    input:
        config_file=ancient(f"results/popsize/{PREFIX}_deminfhelper.yml"),
        vcf=VCF,
    output:
        msmc2_output=f"results/popsize/output_msmc2/{PREFIX}_msmc2.final.txt",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * resources['msmc2']['mem_mb']
    threads:
        resources['msmc2']['threads']
    conda:
        "envs/deminfhelper.yml"
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --msmc2 ;"
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --plot_msmc2"


rule statistics:
    input:
        config_file=ancient(f"results/popsize/{PREFIX}_deminfhelper.yml"),
        vcf=VCF,
    output:
        pca_output=f"results/popsize/output_stats/{PREFIX}_PCA.{PARAMS['plot_format']}",
    conda:
        "envs/deminfhelper.yml"
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --pca ;"
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --gq_distrib"
