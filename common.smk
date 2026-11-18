import re
import gzip

import yaml

deminfhelper_config = "workflow/modules/popsize/config/deminfhelper_config.yml"


def resolve_popsize_params():
    """All module parameters, read from config/deminfhelper_config.yml."""
    with open(deminfhelper_config) as f:
        return yaml.safe_load(f)


def get_popsize_targets(params):
    """(populations, multi_pop, tools) from the resolved params, needed at parse time for rule paths.

    `populations` (name -> member sample IDs) is the single source of truth for
    population(s), same key deminfhelper.py itself reads when run standalone. One
    entry means a single population, named by that entry's key -- byte-identical to
    the module's original single-population behavior. More than one entry splits the
    cohort into that many independently-computed populations, each its own
    results/popsize/<name>/ branch.
    """
    populations = list(params["populations"].keys())
    multi_pop = len(populations) > 1
    tools = {t.strip() for t in params["popsize_tools"].split(",") if t.strip()}
    return populations, multi_pop, tools


def get_population_samples(params, population):
    """Member sample IDs for a population, or None to use every VCF sample ('all')."""
    if population not in params["populations"]:
        raise ValueError(
            f"Population '{population}' is not declared under 'populations' in "
            f"deminfhelper_config.yml."
        )
    ids = params["populations"][population]
    return None if (not ids or ids == "all") else list(ids)


def get_snparcher_paths(config):
    """vcf/ref paths deduced from the snpArcher run, not user-configurable."""
    if config["modules"]["postprocess"]["enabled"]:
        vcf = "results/postprocess/clean_snps.vcf.gz"
    else:
        vcf = "results/vcfs/filtered.vcf.gz"
    return vcf, f"results/reference/{config['reference']['name']}.fa.gz"


def build_deminfhelper_config(out_dir, vcf_file, ref_fasta, population, output, params, samples_override=None):
    # parse VCF header: total genotyped length (from ##contig lines) and sample list
    Nseq = 0
    with gzip.open(vcf_file, 'rt') as vcf:
        line = vcf.readline()
        while line != "" and line.startswith("#"):
            if line[0:8] == "##contig":
                contig_length = int(re.split('[=,]', line)[-1][:-2])
                Nseq += contig_length
            if line.startswith("#CHROM"):
                vcf_samples = line.split()[9:]
                if samples_override is not None:
                    missing = set(samples_override) - set(vcf_samples)
                    if missing:
                        raise ValueError(
                            f"Population sample(s) not found in {vcf_file}: {sorted(missing)}"
                        )
                    samples = samples_override
                else:
                    samples = vcf_samples
                break
            line = vcf.readline()

    samples_str = ", ".join(samples)

    lines = [
        f"out_dir: {out_dir}/",
        "## DATA",
        f"vcf: {vcf_file}",
        f"ref_genome: {ref_fasta}",
        "## L: effective sequence length genotyped",
        f"L: {Nseq}",
        "## Population",
        "populations:",
        f"  {population}: [{samples_str}]",
        "# PARAM",
        f"gen_time: {params['gen_time']}",
        f"mut_rate: {params['mut_rate']}",
        "## PLOTTING",
        f"plot_format: {params['plot_format']}",
        "## COMPUTE SFS",
        f"out_dir_sfs: {out_dir}/",
        f"folded: {params['folded']}",
        f"length_cutoff: {params['length_cutoff']}",
        "## Filtering contigs by name",
        f"contig_filter: {params['contig_filter']}",
        "## SFS if previously computed",
        f"path_to_sfs: {out_dir}/SFS_{population}.fs",
        f"out_dir_stairwayplot2: {out_dir}/output_stairwayplot2/",
        f"summary_file_stw: {out_dir}/output_stairwayplot2/{population}/{population}.final.summary",
        "## Dadi",
        f"out_dir_dadi: {out_dir}/output_dadi/",
        f"optimizations: {params['optimizations']}",
        f"p0: {params['p0']}",
        f"lower_bound: {params['lower_bound']}",
        f"upper_bound: {params['upper_bound']}",
        "## MSMC2",
        f"out_dir_msmc2: {out_dir}/output_msmc2/",
        f"msmc2_kwargs: {params['msmc2_kwargs']}",
        "## SMC++",
        f"out_dir_smcpp: {out_dir}/output_smcpp/",
        f"plot_file_smcpp: {out_dir}/output_smcpp/{population}_inference.csv",
        "## PSMC",
        f"out_dir_psmc: {out_dir}/output_psmc/",
        f"psmc_kwargs: {params['psmc_kwargs']}",
        f"plot_psmc_kwargs: {params['plot_psmc_kwargs']}",
        "## GQ distribution",
        f"out_dir_gq_distrib: {out_dir}/output_stats/",
        "## FINAL INFERENCES",
        f"final_out_dir: {out_dir}/inferences/",
        "## Stats",
        f"out_dir_stats: {out_dir}/output_stats/",
        "## Stairway Plot 2",
        f"path_to_stairwayplot2: {params['path_to_stairwayplot2']}",
        f"blueprint_template: {params['blueprint_template']}",
        "## PCA",
        f"mem: {params['mem']}",
        f"n_clust_kmeans: {params['n_clust_kmeans']}",
    ]

    with open(output, 'w') as config_out:
        config_out.write("\n".join(lines) + "\n")
