#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
inferences.py - Population Demographic Inference Module of DemInfHelper.

This script facilitates demographic inference for populations using various tools and methods.
It provides functions to perform demographic analyses, including estimation of population size histories
and other relevant parameters.

Usage:
1. Prepare input data:
   - Provide VCF (Variant Call Format) files containing genetic data for the population.
   - Define population-specific parameters, such as mutation rate (mu) and generation time (gen_time).

2. Import and use the functions in this script to conduct demographic inference:
   - `run_dadi_cli`: Runs demographic inference using DADI and saves results.
   - `run_stairwayplot2`: Prepares input files for and runs StairwayPlot2 to estimate population size histories.
   - `msmc2`: Performs demographic inference using MSMC2 on multiple contigs.
   - `psmc`: Estimates demographic history using PSMC and generates plots.
   - `smcpp`: Conducts demographic inference using SMC++ on multiple contigs.
   - `run_vcf2smc`: Converts VCF data to SMC++ input format.
   - Other utility functions for data preprocessing and analysis.

3. Customize parameters:
   - Modify parameters within each function as needed, such as population-specific options.

4. Run the script with appropriate arguments to execute demographic inference.

5. Examine the generated output files and plots for insights into the population's demographic history.

Dependencies:
- SMC++
- MSMC2
- PSMC
- DADI
- StairwayPlot2

Note:
- Ensure that the required software tools are installed and accessible.
- Properly configure input data paths and parameters before running the script.
"""

import re
import os
from textwrap import wrap
import subprocess
import itertools
import multiprocessing
import dadi
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import gzip 
import parsing

from numpy import array
from matplotlib import pyplot as plt

def run_dadi_cli(popid, out_dir, sfs_path, optimizations=None,
                 p0 = [1, 1, 10, 5],
                 lower_bounds = [0.1, 0.1, 0.1, 0.1],
                 upper_bounds = [50, 5, 30, 10]):
    """
    Executes dadi-cli for demographic model inference using SFS.
    
    This function runs dadi-cli commands to infer demographic models and
    determine the best fit with the three_epoch model. It allows specifying
    initial parameters, bounds, and optimization steps.
    
    Parameters:
    - popid (str): Identifier for the population or sample.
    - out_dir (str): Output directory path for results.
    - sfs_path (str): Path to the SFS file for model inference.
    - optimizations (int, optional): Number of optimization iterations. If None,
      runs until convergence.
    - p0 (list, optional): Initial model parameter values. Default is [1, 1, 10, 5].
    - lower_bounds (list, optional): Lower bounds for parameters. Default is
      [0.1, 0.1, 0.1, 0.1].
    - upper_bounds (list, optional): Upper bounds for parameters. Default is
      [50, 5, 30, 10].

    The function runs `InferDM` for model inference and then `BestFit` to find
    the best fitting model. Convergence and optimizations are managed based on
    the `optimizations` parameter.

    Outputs:
    - Console output of dadi-cli commands.
    - Result files in the specified output directory.
    """
    if optimizations == None:
        check_convervence_step = 10
    else:
        check_convervence_step = int(optimizations) // 10
    # create dadi file
    if optimizations == None:
        # if no number of optimizations is set, run until convergence.
        cmd1 = "".join(["dadi-cli InferDM --fs ", sfs_path, \
                        " --model three_epoch --lbounds ", " ".join(map(str, lower_bounds)), \
                        " --ubounds ", " ".join(map(str, upper_bounds)),
                        " --p0 ", p0, " --output ", \
                        out_dir+str(popid), " --force-convergence ", str(check_convervence_step) ])
    else:
        cmd1 = "".join(["dadi-cli InferDM --fs ", sfs_path, \
                        " --model three_epoch --lbounds ", " ".join(map(str, lower_bounds)), \
                        " --ubounds ", " ".join(map(str, upper_bounds)), " --output ", \
                        out_dir+str(popid), " --p0 ", " ".join(map(str, p0)), \
                        " --optimizations ", str(optimizations), \
                        # check for convergence every 10 optimizations
                        " --check-convergence ", str(check_convervence_step) ])
    print(cmd1)
    os.system(cmd1)
    cmd2 = "".join(["dadi-cli BestFit --input-prefix ", out_dir+str(popid)+".InferDM", \
                        " --lbounds ", " ".join(map(str, lower_bounds)), \
                        " --ubounds ", " ".join(map(str, upper_bounds)) ])
    print(cmd2)
    os.system(cmd2)    

def input_stairwayplot2(popid, nseq, L, whether_folded, SFS, mu, year_per_generation, stairway_plot_dir, output_path, temp_blueprint):
    """
    Writes the input file for running Stairway Plot v2.

    This function prepares a blueprint file for Stairway Plot v2, which is used 
    for inferring population size history from SFS data.

    Parameters:
    - popid (str): Name of the population.
    - nseq (int): Number of sequences (2*n for diploids).
    - L (int): Total number of sites.
    - whether_folded (bool): Indicates if the VCF is phased (False) or not (True).
    - SFS (list): Site frequency spectrum without monomorphic sites.
    - mu (float): Mutation rate.
    - year_per_generation (float): Generation time in years.
    - stairway_plot_dir (str): Path to the Stairway Plot v2 directory.
    - output_path (str): Path where the output of Stairway Plot v2 will be saved.
    - temp_blueprint (str): Path to the temporary blueprint file.

    The function reads the temporary blueprint file, replaces placeholders with 
    actual values based on the provided parameters, and writes the final 
    blueprint file to the specified output path.
    """

    locals()['project_dir'] = output_path+popid #where the output of stairwayplot2 will be
    locals()['plot_title'] = popid #name of the plots output by stairwayplot2
    with open(temp_blueprint, "r") as temp, open(output_path+str(popid)+".blueprint","w") as out_file:
        line = temp.readline()
        while line != '':
            if line.split(':')[0] in locals().keys():
                if line.split(':')[0] == 'SFS':
                    out_file.write('sfs: ' + ' '.join(map(str, SFS)) + '\n')
                else:
                    out_file.write(line.split(':')[0]+': '+ str(locals().get(line.split(':')[0]))+'\n')
            elif line.split(':')[0] == "nrand" :
                out_file.write('nrand: '+str(int((nseq-2)/4))+' '+ str(int((nseq-2)/2))+' '+ str(int((nseq-2)*3/4))+' '+ str(int(nseq-2))+'\n')
            else:
                out_file.write(line)
            line = temp.readline()

def run_stairwayplot2(popid, out_dir, path_to_stairwayplot2):
    """
    Executes Stairway Plot v2 for population size history inference.

    This function runs Stairway Plot v2 using a blueprint file generated for
    the specified population. It first uses Java to build the stairway plot
    and then executes the generated bash script to complete the process.

    Parameters:
    - popid (str): Identifier for the population.
    - out_dir (str): Output directory path where results and scripts will be saved.
    - path_to_stairwayplot2 (str): Path to the Stairway Plot v2 installation directory.

    The function performs two main steps:
    1. Executes a Java command to build the stairway plot using the blueprint file.
    2. Runs the generated bash script to complete the Stairway Plot v2 analysis.

    Outputs:
    - Console output of the Java command execution.
    - Results of the Stairway Plot v2 analysis are saved in the specified output directory.
    """
    cmd1 = "".join(["java -cp ", path_to_stairwayplot2, " Stairbuilder ", out_dir, popid, ".blueprint"])
    output = subprocess.check_output([cmd1], shell=True)
    output = output.decode(encoding="utf-8")
    print(output)
    cmd2 = "".join([ "bash ", out_dir, popid, ".blueprint.sh"])
    os.system(cmd2)

def run_msmc2_process(args):
    contig, individual, vcf, out_dir, deminfhelper_directory = args
    individual_vcf = f"{out_dir}{contig}_{individual}.vcf.gz"
    # construct and execute the command as in the original function
    cmd = " ".join(["bcftools view -s", individual, "-t", contig, vcf,
                    "| bcftools query - -f '%INFO/DP\n' | awk '{ sum += $1 } END { print sum/NR }' |",
                    "{ read value && minDP=$(echo \"$value / 2\" | bc) && maxDP=$(echo \"$value * 2\" | bc);};",
                    "bcftools view -g ^miss -t", contig, "-s", individual, vcf,
                    "| vcftools --vcf - --minDP $minDP --maxDP $maxDP --recode --stdout | gzip -c >",
                    individual_vcf, ";",
                    "python3", deminfhelper_directory+"/scripts/generate_multihetsep.py", individual_vcf,
                    ">", out_dir+contig+"_"+individual+"_msmc_input.txt"])
    subprocess.run(cmd, shell=True, check=True)
    
def msmc2(contigs, popid, pop_ind, vcf, out_dir, mu, gen_time, num_cpus=None):
    """
    Runs MSMC2 analysis on given contigs with specified parameters.

    This function processes VCF files for each contig to create inputs suitable 
    for MSMC2 analysis. It performs filtering, data conversion, and runs MSMC2.

    Parameters:
    - contigs (list): List of contigs to be analyzed.
    - popid (str): Identifier for the population.
    - pop_ind (int): Index of the population in the VCF file.
    - vcf (str): Path to the VCF file.
    - out_dir (str): Output directory for MSMC2 results.
    - mu (float): Mutation rate.
    - gen_time (float): Generation time.
    - num_cpus (int, optional): Number of CPUs to use. Default is all available.

    The function:
    1. Filters VCF files for each contig based on depth.
    2. Converts filtered VCF files to MSMC2 input format.
    3. Runs MSMC2 on the prepared input files.

    If no contigs are provided or all input files are empty, it raises an error.

    Note:
    - Uses multiprocessing for parallel processing of contigs.
    - Requires bcftools, vcftools, and MSMC2 to be installed and accessible.
    """
    if len(contigs) == 0:
        raise ValueError("Error! No contigs to use! Make sure the threshold matches your data.")
    
    deminfhelper_directory = os.path.dirname(os.path.abspath(__file__))
    num_cpus = num_cpus or multiprocessing.cpu_count()

    pool_args = []
    for contig in contigs:
        for individual in pop_ind:
            pool_args.append((contig, individual, vcf, out_dir, deminfhelper_directory))

    with multiprocessing.Pool(num_cpus) as pool:
        pool.map(run_msmc2_process, pool_args)

    #List only non-empty files for msmc2
    kept_files = [f for f in os.listdir(out_dir) if os.path.isfile(os.path.join(out_dir, f)) and f.endswith("_msmc_input.txt") and os.path.getsize(os.path.join(out_dir, f)) > 0]

    if len(kept_files) == 0:
        raise ValueError("Error! There are no usable file for msmc2, all inputs are empty!")
    # Process each non-empty file
    for filename in kept_files:
        with open(os.path.join(out_dir, filename), 'r') as input_file, open(os.path.join(out_dir, filename + "_temp"), 'w') as output_file:
            for line in input_file:
                # Check if matches the pattern checked by msmc2
                pattern = r"^[A-Za-z0-9_.]+\s\d+\s\d+\s[ACTG01\?,]+$"
                replace_pattern = r"[^A-Za-z0-9_.\sACTG01\?,]"
                # count nb of diff. alleles
                nb_polyall = len(set(line.strip().split("\t")[-1]))
                if nb_polyall == 1:
                    # skip monomorphic sites
                    continue
                if len(line.strip().split("\t")[-1]) > 2:
                    # in MSMC format, the 4th col represents all the alleles of the different samples
                    # as we work sample per sample and chr per chr, we fix this at 2, otherwise, possible inconsistency will make things crash
                    continue
                if not re.match(pattern, line):
                    # Replace hyphens with underscores in the first field
                    if "*" in line:
                        line = line.replace("*", "")
                    line = re.sub(replace_pattern, "_", line)
                    output_file.write(line)
                else:
                    # if line is correct write it as is
                    output_file.write(line)

        # Replace the original file with the corrected one
        os.replace(os.path.join(out_dir, filename + "_temp"), os.path.join(out_dir, filename))
        # Print information or run MSMC2 if needed
        print(f"Processed {filename}.")

    #List only non-empty files for msmc2 after filtering  
    kept_files = [f for f in os.listdir(out_dir) if os.path.isfile(os.path.join(out_dir, f)) and f.endswith("_msmc_input.txt") and os.path.getsize(os.path.join(out_dir, f)) > 0]

    if len(kept_files) == 0:
        raise ValueError("Error! There are no usable file for msmc2, all inputs are empty!")

    cmd5 = " ".join(["msmc2_Linux", ' '.join([os.path.join(out_dir, f) for f in kept_files]), "-o", out_dir+popid+"_msmc2"])
    print(cmd5)
    with open(out_dir+"/"+popid+"_msmc2.log", 'w') as log:
        subprocess.run(cmd5, stdout=log, shell=True)           

def psmc(ref_genome, contigs, popid, pop_ind, vcf, out_dir, mu, gen_time, p="10+22*2+4+6", num_cpus=4):
    """
    Executes the PSMC analysis for population size history inference.

    This function processes VCF files along with a reference genome to generate
    inputs for PSMC and then performs PSMC analysis for each sample in the
    population index.

    Parameters:
    - ref_genome (str): Path to the reference genome. If None, a fake reference
      genome is generated for bcftools consensus.
    - contigs (list): List of contigs to be analyzed.
    - popid (str): Identifier for the population.
    - pop_ind (list): List of indices for the population in the VCF file.
    - vcf (str): Path to the VCF file.
    - out_dir (str): Output directory for PSMC results.
    - mu (float): Mutation rate.
    - gen_time (float): Generation time.
    - p (str, optional): PSMC model parameter string. Default is "10+22*2+4+6",
      specifying 64 atomic time intervals and 28 free interval parameters.
      These parameters are chosen to ensure sufficient recombination events
      are inferred in each interval, avoiding overfitting. You would need to adapt
      these depending on the species you are studying.
    - num_cpus (int, optional): Number of CPUs to use. Default is 4.

    The function processes each contig, generates consensus sequences, and runs
    PSMC, combining results for population-level analysis.

    Outputs:
    - PSMC analysis files in the specified output directory.

    Note:
    - Requires bcftools, seqtk, gzip, and PSMC to be installed and accessible.
    """
    # Get the available CPUs
    available_cpus = list(range(multiprocessing.cpu_count()))
    # Choose the first 'num_cpus' CPUs
    cpu_affinity = available_cpus[:num_cpus]

    if not ref_genome:
        # WIP
        ## will create a fake ref genome for bcftools consensus
        print("Warning: No Reference genome provided. A fake ref genome is being generated for bcftools consensus.")
        
        ref_genome = out_dir+popid+".fa"

        # Initialize variables
        current_chrom = None
        current_pos = 0
        old_pos = 0

        with gzip.open(vcf, 'rt') as vcf_file, open(ref_genome, 'w') as fasta_file:
            line = vcf_file.readline()
            pbar = tqdm(total=0, dynamic_ncols=True, unit='line', unit_scale=True)
            # start parsing the file
            while line:
                if line.startswith("#"):
                    line = vcf_file.readline()
                    continue
                fields = line.split()
                chrom = fields[0]
                pos = int(fields[1])
                if pos == old_pos:
                    line = vcf_file.readline()
                    continue
                ref_allele = fields[3]
                alt_allele = fields[4]
                if chrom != current_chrom:
                    # write header
                    fasta_file.write(f">{chrom}\n")
                    current_chrom = chrom
                    # a CHRM starts at pos 1
                    writing_head = 1
                    # handle indels
                if pos < writing_head:
                    ref_allele = ref_allele[writing_head-pos:]
                    pos = writing_head
                current_pos = pos
                # write N until the current Ref
                for i in range(writing_head, current_pos):
                    fasta_file.write("N")
                    # print(writing_head, current_pos, current_pos-1, ref_allele) 
                    if i % 70 == 0:
                        fasta_file.write("\n")
                fasta_file.write(f"{ref_allele}")
                if (current_pos+len(ref_allele)-1) % 70 == 0:
                    fasta_file.write("\n")
                writing_head = current_pos+len(ref_allele)
                old_pos = pos
                line = vcf_file.readline()
                pbar.update(1)
            pbar.close()

    processes = []
    for sample in pop_ind:
        cmd2 = " ".join(["bcftools consensus -I", vcf, "-s", sample, "-f", ref_genome, "|",
        # read from stdin
        "seqtk seq -F '#'", "-", "|",
        "bgzip >", out_dir+"consensus_"+sample+".fq.gz", ";",
        "fq2psmcfa -q1", out_dir+"consensus_"+sample+".fq.gz", ">", out_dir+sample+"_diploid.psmcfa", ";"
        "psmc -N30 -t15 -r5 -p '"+p+"' -o", out_dir+sample+".psmc", out_dir+sample+"_diploid.psmcfa"])
        print(cmd2)
        with open(out_dir+sample+"_consensus.log", 'w') as log:
            process = subprocess.Popen(cmd2, shell=True, stdout=log)
            # Set CPU affinity
            process_pid = process.pid
            os.sched_setaffinity(process_pid, cpu_affinity)
            processes.append(process)
    for process in processes:
        process.wait()
    cmd3 = " ".join(["cat", out_dir+"*.psmc >", out_dir+'/'+popid+"_combined.psmc.final", ";"
    "psmc_plot.pl", "-g", gen_time, "-x", "0", "-u", mu, "-M '"+",".join(pop_ind)+"'", popid, out_dir+'/'+popid+"_combined.psmc.final" ,";",
    "mv", popid+".*.par", out_dir, "; mv", popid+"*.eps", out_dir])
    print(cmd3)
    with open(out_dir+popid+"_psmc_combine.log", 'w') as log:
        p=subprocess.Popen(cmd3,stdout=log, shell=True)
    p.wait()

def run_vcf2smc(contig, args):
    """
    Converts VCF data to SMC++ input format for a specific contig and population.
    This is a worker meant to be used for multiprocessing in smcpp().
    
    This function takes VCF data and converts it to the SMC++ input format for a
    specified contig and population. It uses the SMC++ 'vcf2smc' command to
    perform this conversion.

    Parameters:
    - contig (str): Name of the contig being processed.
    - args (tuple): A tuple containing the following elements:
        - vcf (str): Path to the VCF file.
        - out_dir (str): Output directory for SMC++ results.
        - popid (str): Identifier for the population.
        - pop_ind (list): List of indices for the population in the VCF file.

    Outputs:
    - SMC++ input file in the specified output directory.
    - Log file with conversion details.

    Note:
    - Requires SMC++ to be installed and accessible.
    """
    vcf, out_dir, popid, pop_ind = args
    cmd1 = ["smc++", "vcf2smc", vcf, out_dir + popid + "_" + contig + ".smc.gz", contig, popid + ":" + ",".join(pop_ind)]
    log_file = out_dir + popid + "_" + contig + "_vcf2smc.log"
    with open(log_file, 'w') as log:
        subprocess.run(cmd1, stdout=log)

def smcpp(contigs, popid, pop_ind, vcf, out_dir, mu, gen_time, num_cpus=None):
    """
    Perform demographic inference using SMC++ on multiple contigs of a population's VCF data.

    This function conducts demographic inference using the SMC++ tool on multiple contigs
    of a population's VCF (Variant Call Format) data. It parallelizes the VCF to SMC++
    conversion process for each contig and then estimates demographic parameters.

    Parameters:
    - contigs (list): List of contig names to process.
    - popid (str): Identifier for the population.
    - pop_ind (list): List of indices for the population in the VCF file.
    - vcf (str): Path to the VCF file containing genetic data.
    - out_dir (str): Output directory for storing SMC++ results and plots.
    - mu (str): Mutation rate for demographic inference.
    - gen_time (str): Generation time for demographic inference.
    - num_cpus (int, optional): Number of CPU cores to use for parallel processing. Default is None.

    Outputs:
    - SMC++ output files and plots in the specified output directory.
    """
    if len(contigs) == 0:
        with open(out_dir + popid + "_model.final.json", 'w') as output:
            output.write("There was an error with SMC++. Please check the logs.")
        raise ValueError("Error! No contigs to use! Make sure the threshold matches your data.")

    # Get the number of CPUs to use
    if num_cpus is None:
        num_cpus = multiprocessing.cpu_count()

    # Parallel processing of contigs
    pool_args = (vcf, out_dir, popid, pop_ind)
    with multiprocessing.Pool(num_cpus) as pool:
        pool.starmap(run_vcf2smc, [(contig, pool_args) for contig in contigs])

    # Once all contig processing is done, proceed with the next steps
    kept_files = [f for f in os.listdir(out_dir) if os.path.isfile(os.path.join(out_dir, f)) and f.endswith(".smc.gz") and os.path.getsize(os.path.join(out_dir, f)) > 0]
    cmd2 = ["smc++", "estimate", "-o", out_dir, "--base", popid, mu] + [os.path.join(out_dir, f) for f in kept_files]
    print(cmd2)
    output2 = subprocess.check_output(cmd2)
    output2 = output2.decode(encoding="utf-8")
    print(output2)

    cmd3 = ["smc++", "plot", out_dir + popid + "_inference.png", out_dir + "/" + popid + ".final.json", "-g", str(gen_time), "-c"]
    subprocess.run(cmd3)



