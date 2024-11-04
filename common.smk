import gzip
from yaml import safe_load

#with open(config["resource_config"], "r") as f:
#    resources = safe_load(f)
#print(resources)

def build_deminfhelper_config(results_folder, vcf_file, prefix, output, global_config, module_config):
    # parse template

    #init vars
    Nseq = 0
    length_cutoff = global_config['contig_size']
    # parse VCF to deduce L and kept samples at this stage
    with gzip.open(vcf_file, 'rt') as vcf:
        line = vcf.readline()
        # parsing the vcf
        while line != "" and line.startswith("#"):
            if line[0:8] == "##contig":
                # keep only contigs that are longer than the length_cutoff parameter
                contig_length = int(re.split('[=,]', line)[-1][:-2])
                #if contig_length >= length_cutoff:
                Nseq += contig_length

            if line.startswith("#CHROM"):
                samples = ','.join(line.split()[9:])
                break
            line = vcf.readline()

    # create config
    with open(module_config, 'r') as module_conf:
        config_parsed = module_conf.read()
    with open(output, 'w') as config_out:
        config_out.write("out_dir: "+results_folder+"/popsize/\n"+\
                         "## DATA\nvcf: "+vcf_file+"\n"+\
                         "ref_genome: "+results_folder+"/data/genome/"+results_folder.split("/")[-1]+".fna\n"\
                         "## L: effective sequence length genotyped\nL: "+str(Nseq)+"\n"+\
                         "## Population(s)\nname_pop: "+prefix+"\n"+\
                         "npop: 1\n# POP 1\n"+prefix+": "+samples+"\n"+\
                         "# PARAM\ngen_time: "+str(global_config['gen_time'])+"\n"+\
                         "mut_rate: "+str(global_config['mut_rate'])+"\n"+\
                         "## COMPUTE SFS\nout_dir_sfs: "+results_folder+"/popsize/\n"+\
                         "folded: "+str(global_config['folded'])+"\n"+\
                         "length_cutoff: "+str(global_config['contig_size'])+"\n"+\
                         "## SFS if previously computed\npath_to_sfs: "+results_folder+"/popsize/SFS_"+prefix+".fs\n"+\
                         "out_dir_stairwayplot2: "+results_folder+"/popsize/output_stairwayplot2/\n"+\
                         "summary_file_stw: "+results_folder+"/popsize/output_stairwayplot2/"+prefix+"/"+prefix+".final.summary\n"+\
                         "out_dir_dadi: "+results_folder+"/popsize/output_dadi/\n"+\
                         "optimizations: "+str(global_config['dadi_optimizations'])+"\n"+\
                         "out_dir_msmc2: "+results_folder+"/popsize/output_msmc2/\n"+\
                         "## SMC++\nout_dir_smcpp: "+results_folder+"/popsize/output_smcpp/\n"+\
                         "## PSMC\nout_dir_psmc: "+results_folder+"/popsize/output_psmc/\n"+\
                         "plot_file_smcpp: "+results_folder+"/popsize/output_smcpp/"+prefix+"_inference.csv\n"+\
                         "## GQ distribution\nout_dir_gq_distrib: "+results_folder+"/popsize/output_stats/\n"+\
                         "## FINAL INFERENCES\nfinal_out_dir: "+results_folder+"/popsize/inferences/\n"+ \
                         "## Stats\nout_dir_stats: "+results_folder+"/popsize/output_stats/\n" + \
                         "## PCA K-means clustering\nn_clust_kmeans: 1\n")
        config_out.write(config_parsed)
