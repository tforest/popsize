import gzip

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
                        "## L: effective sequence length genotyped\nL: "+str(Nseq)+"\n"+\
                        "## Population(s)\nname_pop: "+prefix+"\n"+\
                        "npop: 1\n# POP 1\n"+prefix+": "+samples+"\n"+\
                        "# PARAM\ngen_time: "+str(global_config['gen_time'])+"\n"+\
                        "mut_rate: "+str(global_config['mut_rate'])+"\n"+\
                        "## COMPUTE SFS\nout_dir_sfs: "+results_folder+"/popsize/\n"+\
                        "folded: "+str(global_config['folded'])+"\n"+\
                        "## SFS if previously computed\npath_to_sfs: "+results_folder+"/popsize/SFS_"+prefix+".txt\n"+\
                        "out_dir_stairwayplot2: "+results_folder+"/popsize/output_stairwayplot2/\n"+\
                        "summary_file_stw: "+results_folder+"/output_stairwayplot2/"+prefix+"/"+prefix+".final.summary\n"+\
                        "## SMC++\nout_dir_smcpp: "+results_folder+"/popsize/output_smcpp/\n"+\
                        "length_cutoff: 100000\n"+\
                        "plot_file_smcpp: "+results_folder+"/popsize/output_smcpp/"+prefix+"_inference.csv\n"+\
                        "## GQ distribution\nout_dir_gq_distrib: "+results_folder+"/output_gq_distrib/\n"+\
                        "## FINAL INFERENCES\nfinal_out_dir: "+results_folder+"/inferences/\n")
        config_out.write(config_parsed)
