include: "common.smk"
include: "config/resources.smk"

PARAMS = resolve_popsize_params()
POPULATIONS, MULTI_POP, _TOOLS = get_popsize_targets(PARAMS)
VCF, REF_FASTA = get_snparcher_paths(config)

# Bare wildcard tokens for f-string embedding. _POP is always a real "{population}"
# wildcard (even with a single population, trivially constrained to one legal value,
# so behavior/paths are unaffected). _POP_DIR additionally nests output under a
# per-population directory, but only when there truly is more than one population --
# this keeps the single-population file layout byte-identical to before.
_POP = "{population}"
_POP_DIR = f"{_POP}/" if MULTI_POP else ""

wildcard_constraints:
    population = "|".join(re.escape(p) for p in POPULATIONS)

_TARGET_PATTERNS = [p for p in [
    f"results/popsize/{_POP_DIR}SFS_{_POP}.fs",
    f"results/popsize/{_POP_DIR}output_stats/{_POP}_PCA.{PARAMS['plot_format']}",
    f"results/popsize/{_POP_DIR}output_smcpp/{_POP}.final.json" if "smcpp" in _TOOLS else None,
    f"results/popsize/{_POP_DIR}output_stairwayplot2/{_POP}/{_POP}.final.summary" if "swp2" in _TOOLS else None,
    f"results/popsize/{_POP_DIR}output_dadi/{_POP}.InferDM.bestfits" if "dadi" in _TOOLS else None,
    f"results/popsize/{_POP_DIR}output_psmc/{_POP}.eps" if "psmc" in _TOOLS else None,
    f"results/popsize/{_POP_DIR}output_msmc2/{_POP}_msmc2.final.txt" if "msmc2" in _TOOLS else None,
    f"results/popsize/{_POP_DIR}{_POP}_combined_plot.{PARAMS['plot_format']}",
] if p is not None]

_ALL_TARGETS = [f for pat in _TARGET_PATTERNS for f in expand(pat, population=POPULATIONS)]

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
        config_file=f"results/popsize/{_POP_DIR}{_POP}_deminfhelper.yml",
    run:
        out_dir = f"results/popsize/{wildcards.population}" if MULTI_POP else "results/popsize"
        build_deminfhelper_config(
            out_dir=out_dir,
            vcf_file=input.vcf,
            ref_fasta="results/popsize/ref.fa",
            population=wildcards.population,
            output=output.config_file,
            params=PARAMS,
            samples_override=get_population_samples(PARAMS, wildcards.population),
        )


rule compute_sfs:
    input:
        config_file=ancient(f"results/popsize/{_POP_DIR}{_POP}_deminfhelper.yml"),
    output:
        sfs_file=f"results/popsize/{_POP_DIR}SFS_{_POP}.fs",
    conda:
        "envs/deminfhelper.yml"
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --sfs --config_file {input.config_file} ;"
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --plot_sfs"


rule stairwayplot2:
    input:
        config_file=ancient(f"results/popsize/{_POP_DIR}{_POP}_deminfhelper.yml"),
        sfs_file=f"results/popsize/{_POP_DIR}SFS_{_POP}.fs",
    output:
        stairwayplot2_summary=f"results/popsize/{_POP_DIR}output_stairwayplot2/{_POP}/{_POP}.final.summary",
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
        config_file=ancient(f"results/popsize/{_POP_DIR}{_POP}_deminfhelper.yml"),
        vcf_index=f"{VCF}.tbi",
    output:
        smcpp_summary=f"results/popsize/{_POP_DIR}output_smcpp/{_POP}.final.json",
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
        config_file=ancient(f"results/popsize/{_POP_DIR}{_POP}_deminfhelper.yml"),
        sfs_file=f"results/popsize/{_POP_DIR}SFS_{_POP}.fs",
    output:
        dadi_summary=f"results/popsize/{_POP_DIR}output_dadi/{_POP}.InferDM.bestfits",
    conda:
        "envs/deminfhelper.yml"
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --dadi ;"
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --plot_dadi"


rule psmc:
    input:
        config_file=ancient(f"results/popsize/{_POP_DIR}{_POP}_deminfhelper.yml"),
        vcf=VCF,
        vcf_index=f"{VCF}.tbi",
        ref_genome="results/popsize/ref.fa",
    output:
        psmc_output=f"results/popsize/{_POP_DIR}output_psmc/{_POP}_combined.psmc.final",
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
        config_file=ancient(f"results/popsize/{_POP_DIR}{_POP}_deminfhelper.yml"),
        psmc_output=f"results/popsize/{_POP_DIR}output_psmc/{_POP}_combined.psmc.final",
    output:
        psmc_figure=f"results/popsize/{_POP_DIR}output_psmc/{_POP}.eps",
    params:
        plot_psmc_kwargs=PARAMS["plot_psmc_kwargs"],
    conda:
        "envs/deminfhelper.yml"
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} "
        "--plot_psmc_kwargs '{params.plot_psmc_kwargs}' --plot_psmc"


rule combined_plot:
    input:
        config_file=ancient(f"results/popsize/{_POP_DIR}{_POP}_deminfhelper.yml"),
        dadi=f"results/popsize/{_POP_DIR}output_dadi/{_POP}.InferDM.bestfits" if "dadi" in _TOOLS else [],
        smcpp=f"results/popsize/{_POP_DIR}output_smcpp/{_POP}.final.json" if "smcpp" in _TOOLS else [],
        swp2=f"results/popsize/{_POP_DIR}output_stairwayplot2/{_POP}/{_POP}.final.summary" if "swp2" in _TOOLS else [],
        msmc2=f"results/popsize/{_POP_DIR}output_msmc2/{_POP}_msmc2.final.txt" if "msmc2" in _TOOLS else [],
        psmc=f"results/popsize/{_POP_DIR}output_psmc/{_POP}.eps" if "psmc" in _TOOLS else [],
    output:
        combined_plot=f"results/popsize/{_POP_DIR}{_POP}_combined_plot.{PARAMS['plot_format']}",
    conda:
        "envs/deminfhelper.yml"
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --combined_plot"


rule msmc2:
    input:
        config_file=ancient(f"results/popsize/{_POP_DIR}{_POP}_deminfhelper.yml"),
        vcf=VCF,
    output:
        msmc2_output=f"results/popsize/{_POP_DIR}output_msmc2/{_POP}_msmc2.final.txt",
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
        config_file=ancient(f"results/popsize/{_POP_DIR}{_POP}_deminfhelper.yml"),
        vcf=VCF,
    output:
        pca_output=f"results/popsize/{_POP_DIR}output_stats/{_POP}_PCA.{PARAMS['plot_format']}",
    conda:
        "envs/deminfhelper.yml"
    shell:
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --pca ;"
        "python3 workflow/modules/popsize/scripts/deminfhelper.py --config_file {input.config_file} --gq_distrib"
