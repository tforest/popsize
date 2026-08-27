import re
import gzip

import yaml

_TEMPLATE_PATH = "workflow/modules/popsize/config/deminfhelper_template.yml"


def resolve_popsize_params():
    """All module parameters, read from config/deminfhelper_template.yml."""
    with open(_TEMPLATE_PATH) as f:
        return yaml.safe_load(f)


def get_popsize_targets(params):
    """(prefix, tools) from the resolved params, needed at parse time for rule paths."""
    prefix = params["pop_name"]
    tools = {t.strip() for t in params["popsize_tools"].split(",") if t.strip()}
    return prefix, tools


def get_snparcher_paths(config):
    """vcf/ref paths deduced from the snpArcher run, not user-configurable."""
    return "results/vcfs/raw.vcf.gz", f"results/reference/{config['reference']['name']}.fa.gz"


def build_deminfhelper_config(out_dir, vcf_file, ref_fasta, prefix, output, params):
    # parse VCF header: total genotyped length (from ##contig lines) and sample list
    Nseq = 0
    with gzip.open(vcf_file, 'rt') as vcf:
        line = vcf.readline()
        while line != "" and line.startswith("#"):
            if line[0:8] == "##contig":
                contig_length = int(re.split('[=,]', line)[-1][:-2])
                Nseq += contig_length
            if line.startswith("#CHROM"):
                samples = ','.join(line.split()[9:])
                break
            line = vcf.readline()

    lines = [
        f"out_dir: {out_dir}/",
        "## DATA",
        f"vcf: {vcf_file}",
        f"ref_genome: {ref_fasta}",
        "## L: effective sequence length genotyped",
        f"L: {Nseq}",
        "## Population(s)",
        f"name_pop: {prefix}",
        "npop: 1",
        "# POP 1",
        f"{prefix}: {samples}",
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
        f"path_to_sfs: {out_dir}/SFS_{prefix}.fs",
        f"out_dir_stairwayplot2: {out_dir}/output_stairwayplot2/",
        f"summary_file_stw: {out_dir}/output_stairwayplot2/{prefix}/{prefix}.final.summary",
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
        f"plot_file_smcpp: {out_dir}/output_smcpp/{prefix}_inference.csv",
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
