#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DemInfHelper v0.1.

Authors:
Thomas Forest (thomas.forest1@edu.mnhn.fr)
Swan Portalier
Camille Steux


DemInfHelper is a versatile command-line tool for performing population genetics analysis, including the computation of the Site Frequency Spectrum (SFS) and various demographic inference methods. It is designed to facilitate the analysis of genetic data, helping researchers understand the demographic history of populations.

Usage:
1. Install the required modules and dependencies.
2. Run DemInfHelper with various command-line arguments to perform specific tasks.

## Key Features:
- Compute the SFS from a VCF (Variant Call Format) file.
- Transform and normalize the SFS.
- Perform demographic inference using methods like StairwayPlot2, dadi, MSMC2, PSMC, and SMC++.
- Plot demographic inferences and SFS distributions.
- Analyze genotype quality (GQ) distributions.
- Conduct Principal Component Analysis (PCA) on genetic data.

## Command-Line Arguments:
- `--config_file`: Path to the configuration file (optional).
- `--cpus`: Number of CPU threads to use (optional).
- `--sfs`: Compute the SFS from the VCF file.
- `--sfs_transformed`: Normalize and transform the SFS.
- `--plot_sfs`: Plot the SFS.
- `--stairwayplot2`: Run StairwayPlot2 for demographic inference.
- `--plot_stairwayplot2`: Plot StairwayPlot2 results.
- `--dadi`: Run demographic inference using dadi.
- `--plot_dadi`: Create population size plots from dadi output.
- `--msmc2`: Run MSMC2 for demographic inference.
- `--plot_msmc2`: Plot MSMC2 results.
- `--psmc`: Run PSMC for demographic inference.
- `--plot_psmc`: Plot PSMC inference results.
- `--gq_distrib`: Compute the Genotype Quality (GQ) distribution.
- `--smcpp`: Run SMC++ for demographic inference.
- `--plot_smcpp`: Plot SMC++ results.
- `--Gplot`: Plot all inferences on the same graph.
- `--folded`: Fold the SFS (default: True).
- `--pca`: Compute Principal Component Analysis (PCA).
- `--plot_pca`: Plot PCA results.
- `--n_clust_kmeans`: Number of clusters for k-means in PCA.
- Additional configuration parameters.

## Usage:
- Run DemInfHelper with appropriate arguments to perform desired analysis tasks :

python3 DemInfHelper/deminfhelper.py --config_file config.txt --sfs
python3 DemInfHelper/deminfhelper.py --config_file config.txt --pca --n_clust_kmeans 3
python3 DemInfHelper/deminfhelper.py --config_file config.txt --stairwayplot2
python3 DemInfHelper/deminfhelper.py --config_file config.txt --plot_stairwayplot2

## Directory Structure:
- The tool creates directories for output and results.

"""

## IMPORT MODULES

import gzip
import matplotlib.pyplot as plt
import argparse
import numpy as np
import os


## IMPORT CUSTOM FUNCTIONS
if __package__ is None or __package__ == '':
    from parsing import *
    from inferences import *
    from sfs import *
    from plots import *
else:
    from .parsing import * 
    from .inferences import *
    from .sfs import *
    from .plots import *
    
def parse_args():
    ## CMD ARGUMENTS
    parser = argparse.ArgumentParser(description='Computes the sfs from the vcf and runs demography inference softwares.')
    #mandatory arguments
    parser.add_argument("--config_file", help="path to the configuration file")
    #optional arguments
    parser.add_argument("--cpus", help="# CPU threads to use",  type=int)
    #SFS
    parser.add_argument("--sfs", help = "to compute the sfs", action = "store_true")
    parser.add_argument("--sfs_transformed", help = "to normalize the sfs", action = "store_true")
    parser.add_argument("--plot_sfs", help = "to plot the sfs", action = "store_true")
    #Stairwayplot2
    parser.add_argument("--stairwayplot2", help = "to run stairwayplot2", action = "store_true")
    parser.add_argument("--plot_stairwayplot2", help = "to run stairwayplot2", action = "store_true")
    #Dadi
    parser.add_argument("--dadi", help = "to run dadi: the sfs must not be transformed", action = "store_true")
    parser.add_argument("--plot_dadi", help = "to create popsize plot from dadi output.", action = "store_true")
    #MSMC2
    parser.add_argument("--msmc2", help = "to run msmc2", action = "store_true")
    parser.add_argument("--plot_msmc2", help = "to plot msmc2", action = "store_true")
    #PSMC
    parser.add_argument("--psmc", help = "to run PSMC", action = "store_true")
    parser.add_argument("--plot_psmc", help = "to plot psmc inference", action = "store_true")
    #PL distribution
    parser.add_argument("--gq_distrib", help = "to compute the GQ (genotype quality) distribution", action = "store_true")
    #SMCPP
    parser.add_argument("--smcpp", help = "run smcpp", action = "store_true")
    parser.add_argument("--plot_smcpp", help = "to plot smcpp inference", action = "store_true")
    parser.add_argument("--Gplot", help = "to plot all inferences on the same graph", action = "store_true")
    parser.add_argument("--folded", help = "Fold the SFS. Default: True", action = "store_true", default=True)
    # Statistics
    # PCA
    parser.add_argument("--pca", help = "Compute PCA using Plink2", action = "store_true")
    parser.add_argument("--plot_pca", help = "Compute PCA using Plink2", action = "store_true")
    parser.add_argument("--n_clust_kmeans", help="Defines the number of k clusters for the k-means on the PCA",  type=int)
    # Config args; if no config file
    parser.add_argument("--popid", help="a population identifier; eg. species name",  type=str)
    parser.add_argument("--samples", help="a list of samples to use in the VCF. By default all samples are taken.", default="all", type=str)
    parser.add_argument("--vcf", help="the vcf file path, only gz format",  type=str)
    parser.add_argument("--gentime", help="the generation time of the species in years. Eg 1.5", type=float)
    parser.add_argument("--mu", help="the average mutation rate per site per generation")
    parser.add_argument("--out", help="Output path of the analysis. Will create a dir.",  type=str)
    parser.add_argument("--out_dir_stats", help="Output path of the stats. Will create a dir.",  type=str)
    parser.add_argument("--L", help="The actual size of the genotyped sequence length that produced the VCF, before any filters.",  type=int)
    parser.add_argument("--n", help="Number of sequences that produced the VCF. Eg for 5 dipl individuals, n=10.",  type=int)

    args = parser.parse_args()
    # if args.sfs:
    #     optional_args = [args.popid, args.vcf, args.out]
    # else:
    #     optional_args = [args.popid, args.gentime, args.mu, args.out]

    # if not args.config_file and not all(optional_args):
    #     parser.error("At least one of the optional arguments is missing: popid, samples, vcf, gentime, mu, out, L, n.\nMaybe you forgot to provide a config file?")

    return args

def main():
    # parse args
    args = parse_args()

    ## CONFIG FILE
    if args.popid is None:
        param = parse_config(args.config_file)
    else:
        program_path = "/".join(os.path.abspath(__file__).split("/")[:-1])+"/"
        param = {
            'out_dir': args.out,
            'vcf': args.vcf,
            # for now only use as a single pop.
            'name_pop': [args.popid],
            'npop': 1,
            args.popid: args.samples,
            'folded': args.folded,
            'gen_time': args.gentime,
            'mut_rate': args.mu,
            'out_dir_sfs': args.out+'/output_sfs/',
            'path_to_sfs': args.out+'/output_sfs/SFS_'+args.popid+'.fs',
            'path_to_stairwayplot2': program_path+'/bin/stairway_plot_es/',
            'blueprint_template': program_path+'/bin/template.blueprint',
            'out_dir_stairwayplot2': args.out+'/output_stairwayplot2/',
            'summary_file_stw': args.out+'/output_stairwayplot2/hirrus/hirrus.final.summary',
            'L': args.L,
            # Dadi params
            'lower_bound': '1, 1, 0.05, 0.01',
            'p0': '0.01, 0.001, 0.01, 0.01',
            'upper_bound': '10, 4, 0.1, 10',
            'optimizations': '100',
            'out_dir_dadi': args.out+'/output_dadi/',
            'out_dir_smcpp': args.out+'/output_smcpp/',
            'plot_file_smcpp': args.out+'/output_smcpp/hirrus_inference.csv',
            'out_dir_gq_distrib': args.out+'/output_stats/',
            'out_dir_stats': args.out+'/output_stats/',
            'final_out_dir': args.out+'/inferences/',
            # default length of contig to keep, useful for SMC++
            'length_cutoff': 100000,
            'ref_genome': None,
            'n_clust_kmeans': args.n_clust_kmeans,
            'cpus': 1
        }
        for p in param["name_pop"]:
            param[p] = param[p].split(",")
            param["n_"+p] = len(param[p])
    
    param["cpus"] = args.cpus
    # Add args to the config for hybrid config + parameters set in args
#    for arg_name in vars(args):
#        arg_value = getattr(args, arg_name)
#        if arg_name not in param:
#            # config value is kept over args
#            param[arg_name] = arg_value
    ## CREATING DIRECTORIES
    if not os.path.exists(param["out_dir"]):
        os.makedirs(param["out_dir"])

    if args.smcpp or args.stairwayplot2 or args.dadi:
        if not os.path.exists(param["final_out_dir"]):
            os.makedirs(param["final_out_dir"])

    # Compute the SFS
    if args.sfs or args.gq_distrib:
        res_pars = vcf_line_parsing(PARAM = param, SFS = args.sfs, SMCPP = args.smcpp, GQ = args.gq_distrib)
        if not param['L']:
            # Needed estimated number of sequence genotyped.
            # from GADMA
            # Assume total length of sequence that was used for SFS building is equal to Nseq.
            # From this data total number of X SNP’s were received.
            # But not all of them were used for SFS: some of them (e.g. total number of Y SNP’s) were filtered out.
            # Then we should count filtered SNP’s and take L value the following way:
            # L = (X - Y) / X * Nseq
            param["L"] = res_pars[2]
            print("Computed L=", param["L"])
        if args.config_file:
            # add/update L_computed in the config file
            param["L_computed"] = res_pars[2]
            print("Adding L_computed=",  param["L_computed"], "to", args.config_file)
            update_config(config_dict = param, 
                          config_file = args.config_file)
    if args.sfs:
        if not os.path.exists(param["out_dir_sfs"]):
            # create the sfs output directory if it does not exists
            os.makedirs(param["out_dir_sfs"])
        SFS_dict = res_pars[0]
        for p in param["name_pop"]:
            # Build coverage plot with current filter
            genotyping_coverage_plot(popid = p, snp_coverage = res_pars[3],
                                     out_dir_stats = param["out_dir_stats"],
                                     filter_prefix = None, nb_bins = 10)
            if param["folded"]:
                # if folded is set to True, store in a string the SFS state
                # for later use. 
                folded_string = "folded"
                # Warning! n_bins*2 only for diploids
                n_bins = len(SFS_dict[p]) * 2
            else:
                folded_string = "unfolded"
                n_bins = len(SFS_dict[p])
            with open(param["out_dir_sfs"]+"SFS_"+p+".fs", 'w') as sfs_out:
                sfs_out.write(str(n_bins)+" "+folded_string+" "+'"'+p+'"\n')
                sfs_out.write(" ".join(map(str, SFS_dict[p])))
                if param["folded"]:
                    sfs_out.write(" 0"*len(SFS_dict[p])+"\n")
                else:
                    # add the trailing carriage return
                    sfs_out.write("\n")
                # generate the mask for a valid dadi.Fs spectrum
                for sfs_bin in range(n_bins):
                    if sfs_bin == 0:
                        # Mask the sfs[0] monomorphic sites
                        sfs_out.write("1 ")
                    elif (param["folded"] == False and sfs_bin == len(SFS_dict[p])):
                        # Mask if it the last bin of the SFS and unfolded sfs
                        sfs_out.write("1 ")
                    elif param["folded"] == True and sfs_bin > (n_bins / 2)-1:
                        # Mask all sites that are null in the SFS for folded spectra
                        sfs_out.write("1 ")
                    else:
                        # Part that is not masked
                        sfs_out.write("0 ")
                # trailing carriage return
                sfs_out.write("\n")
                
    if args.sfs_transformed:
        if args.sfs == False:
            SFS_dict = {}
            for p in param["name_pop"]:
                if "path_to_sfs" not in param.keys():
                    print("--sfs flag or path_to_sfs missing")
                else:
                    with open(param["path_to_sfs"], "rt") as sfs:
                        line=sfs.readline()
                        while line != "":
                            SFS_dict[p] = [int(i) for i in line[:-1].split(",")]
                            line = sfs.readline()
        SFS_dict_trans = {}
        for p in param["name_pop"]:
            SFS_dict_trans[p] = transform_sfs(sfs = SFS_dict[p], n = param["n_"+p], \
                    folded = param["folded"])
            with open(param["out_dir_sfs"]+"SFS_"+p+"_transformed.txt", 'w') as sfs_out:
                sfs_out.write(",".join(map(str, SFS_dict_trans[p])))
            if args.plot_sfs:
                plot_sfs(sfs = SFS_dict_trans[p], plot_title = "SFS (transformed) "+p, output_file = param["out_dir_sfs"]+p+"_trans.png")

    # Run Stairwayplot2
    if args.stairwayplot2:
        if not os.path.exists(param["out_dir_stairwayplot2"]):
            os.makedirs(param["out_dir_stairwayplot2"])
        if args.sfs == False:
            if "path_to_sfs" not in param.keys():
                print("--sfs flag or path_to_sfs missing")
            else:
                SFS_dict = {}
                sfs_list = parse_sfs(param["path_to_sfs"])
                SFS_dict[param["name_pop"][0]] = sfs_list

        for p in param["name_pop"]:
            if param["folded"]:
                nseq = len(SFS_dict[p])*2
            else:
                nseq = len(SFS_dict[p])
            input_stairwayplot2(popid = p, nseq = nseq, L= param["L"], whether_folded = param["folded"], \
                            SFS = SFS_dict[p] , mu = param["mut_rate"], year_per_generation = param["gen_time"], \
                            stairway_plot_dir = param["path_to_stairwayplot2"], output_path = param["out_dir_stairwayplot2"], \
                            temp_blueprint = param["blueprint_template"])
            run_stairwayplot2(popid = p, out_dir = param["out_dir_stairwayplot2"], path_to_stairwayplot2 = param["path_to_stairwayplot2"])
        if args.plot_stairwayplot2:
            for p in param["name_pop"]:
                plot_stairwayplot2(popid = p, summary_file = "".join([param["out_dir_stairwayplot2"], p, "/", p,".final.summary"]), \
                               out_dir = param["final_out_dir"])

    # Run dadi
    if args.dadi:
        if not os.path.exists(param["out_dir_dadi"]):
            os.makedirs(param["out_dir_dadi"])
        if args.sfs == False:
            for p in param["name_pop"]:
                if "path_to_sfs" not in param.keys():
                    raise ValueError("--sfs flag or path_to_sfs missing")

        for p in param["name_pop"]:
            run_dadi_cli(popid = p, out_dir = param["out_dir_dadi"],
                         sfs_path = param["out_dir_sfs"]+"SFS_"+p+".fs",
                         p0 = eval(param["p0"]),
                         lower_bounds = eval(param["lower_bound"]),
                         upper_bounds = eval(param["upper_bound"]),
                         optimizations = param["optimizations"])

    # Dadi plot
    if args.plot_dadi:
        for p in param["name_pop"]:
            dadi_vals = dadi_output_parse(param["out_dir_dadi"]+p+".InferDM.bestfits")
            plot_dadi_output_three_epochs(dadi_vals,p,out_dir = param['final_out_dir'],
                                          mu = eval(param['mut_rate']),
                                          gen_time = eval(param['gen_time']),
                                          L = eval(param['L']))
    #GQ distribution
    if args.gq_distrib:
        if not os.path.exists(param["out_dir_gq_distrib"]):
            os.makedirs(param["out_dir_gq_distrib"])
        GQ_dict = res_pars[1]
        for p in param["name_pop"]:
            plot_distrib_gq(popid = p, gq = GQ_dict[p], out_dir_gq = param["out_dir_gq_distrib"] )
    # PCA
    if args.pca:
        if not os.path.exists(param["out_dir_stats"]):
            os.makedirs(param["out_dir_stats"])
        for p in param["name_pop"]:
            pca_from_vcf(popid = p, vcf_file = param["vcf"],
                         nb_samples = param["n_"+p],
                         out_dir = param["out_dir_stats"])
    ##SMC++
    if args.smcpp:
        contigs = get_contigs_lengths(vcf = param["vcf"], length_cutoff=param["length_cutoff"])
        if not os.path.exists(param["out_dir_smcpp"]):
            os.makedirs(param["out_dir_smcpp"])
        for p in param["name_pop"]:
            smcpp(contigs = contigs, popid = p, pop_ind = param[p], vcf = param["vcf"], \
                  out_dir = param["out_dir_smcpp"], mu = param["mut_rate"], 
                  gen_time = param["gen_time"], num_cpus=param["cpus"])
    ##MSMC2
    if args.msmc2:
        contigs = get_contigs_lengths(vcf = param["vcf"], length_cutoff = param["length_cutoff"])
        if not os.path.exists(param["out_dir_msmc2"]):
            os.makedirs(param["out_dir_msmc2"])
        for p in param["name_pop"]:
            msmc2(contigs = contigs, popid = p, pop_ind = param[p], vcf = param["vcf"], \
                  out_dir = param["out_dir_msmc2"], mu = param["mut_rate"], gen_time = param["gen_time"],
                  num_cpus=param["cpus"])
    ##PSMC
    if args.psmc:
        contigs = get_contigs_lengths(vcf = param["vcf"], length_cutoff = param["length_cutoff"])
        if not os.path.exists(param["out_dir_psmc"]):
            os.makedirs(param["out_dir_psmc"])
        for p in param["name_pop"]:
            psmc(ref_genome = param["ref_genome"], contigs = contigs, popid = p, pop_ind = param[p], vcf = param["vcf"], \
                 out_dir = param["out_dir_psmc"], mu = param["mut_rate"], gen_time = param["gen_time"])
    # Plot SFS
    if args.plot_sfs:
        SFS_dict = {}
        for p in param["name_pop"]:
            sfs_list = parse_sfs(param["path_to_sfs"])
            SFS_dict[param["name_pop"][0]] = sfs_list
            plot_sfs(sfs = SFS_dict[p], plot_title = "SFS "+p, output_file = param["out_dir_sfs"]+"SFS_"+p+".png")

    # Plot PCA
    if args.plot_pca:
        for p in param["name_pop"]:
            plot_pca(plink_eigenvec=param["out_dir_stats"]+"/plink/"+p+".pca.eigenvec",
                     plink_eigenval=param["out_dir_stats"]+"/plink/"+p+".pca.eigenval",
                     popid=p, out_dir=param["out_dir_stats"], n_clusters=param["n_clust_kmeans"])

    # Plot StairwayPlot2
    if args.plot_stairwayplot2 and args.stairwayplot2==False:
        if not os.path.exists(param["out_dir_stairwayplot2"]):
            os.makedirs(param["out_dir_stairwayplot2"])
        for p in param["name_pop"]:
            if "".join(["summary_file_stw"]) not in param.keys():
                print("path to the population final summary file missing")
            else:
                plot_stairwayplot2(popid = p, summary_file = param["summary_file_stw"], \
                               out_dir = param["final_out_dir"])
    if args.plot_msmc2:
        for p in param["name_pop"]:
            plot_msmc2(popid = p, summary_file = "".join([param["out_dir_msmc2"], "/", p, "_msmc2.final.txt"]), \
                           mu = param["mut_rate"], gen_time = param["gen_time"], out_dir = param["final_out_dir"])
            
    ## Plot PSMC
    if args.plot_psmc:
        for p in param["name_pop"]:
            plot_psmc(popid = p, sample_names = param[p],
                      psmc_output_file = param["out_dir_psmc"]+"/"+p+"_combined.psmc.final",
                      plot_output_prefix = param["out_dir_psmc"]+"/"+p,
                      gen_time=param["gen_time"], mut_rate=param["mut_rate"],
                      out_dir = param["final_out_dir"])
    # Plot SMC++
    if args.plot_smcpp:
        for p in param["name_pop"]:
            plot_smcpp(popid = p, summary_file = param["plot_file_smcpp"], out_dir = param["final_out_dir"])


    ##Gplot
    if args.Gplot:
        for p in param["name_pop"]:
            Gplot(T_scaled_gen=param["out_dir_dadi"]+"popt_"+p+"_dadi.txt",gen_time=param["gen_time"],dadi_vals_list=dadi_output_parse(param["out_dir_dadi"]+"output_"+p+".dadi"),out_dir=param["out_dir"], title =p,name_pop=p,popid = p, summary_file2 = param["plot_file_smcpp"],summary_file = param["summary_file_stw"])


if __name__ == "__main__":
    main()
