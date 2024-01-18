#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Creates the input files to run the inferences softwares and run the softwares
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
    #writes the input file to run stairwayplot2
    #popid: name of the population
    #nseq: number of sequences (2*n for diploids)
    #L: total number of sites?
    #whether_folded: if the vcf is phased or not
    #SFS: sfs as a list without the monomorphic sites
    #mu: mutation rate
    #year_per_generation: generation time
    #stairway_plot_dir: path to stairwayplot2
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
    cmd1 = "".join(["java -cp ", path_to_stairwayplot2, " Stairbuilder ", out_dir, popid, ".blueprint"])
    output = subprocess.check_output([cmd1], shell=True)
    output = output.decode(encoding="utf-8")
    print(output)
    cmd2 = "".join([ "bash ", out_dir, popid, ".blueprint.sh"])
    os.system(cmd2)

def msmc2(contigs, popid, pop_ind, vcf, out_dir, mu, gen_time, num_cpus=None):
    if len(contigs) == 0:
        print("Error! No contigs to use! Make sure the threshold matches your data.")
    deminfhelper_directory = os.path.dirname(os.path.abspath(__file__))

    # Get the available CPUs
    available_cpus = list(range(multiprocessing.cpu_count()))
    # Choose the first 'num_cpus' CPUs
    cpu_affinity = available_cpus[:num_cpus]
    
    processes = []  # List to store Process objects

    for contig in contigs:
        # Process each contig
        cmd2 = " ".join(["bcftools", "view", "-t", contig, vcf, "|",
                         "bcftools", "query", "-", "-f", "'%INFO/DP\n'", "|",
                         "awk '{ sum += $1 } END { print sum/NR }' | ",
                         '{ read value && minDP=$(echo "$value / 2" | bc) && maxDP=$(echo "$value * 2" | bc);};',
                         "bcftools view -g ^miss -t", contig, vcf,
                         "|vcftools --vcf - --minDP $minDP --maxDP $maxDP",
                         "--recode --stdout | gzip -c >", out_dir+contig+".vcf.gz ;",
                         "python3", deminfhelper_directory+"/scripts/generate_multihetsep.py",
                         out_dir+contig+".vcf.gz", ">", out_dir+contig+"_msmc_input.txt"])
        print(cmd2)
        p = multiprocessing.Process(target=subprocess.run, args=(cmd2,), kwargs={'shell': True, 'check': True})
        processes.append(p)
        p.start()
    # Wait for all processes to finish
    for process in processes:
        process.join()

    #List only non-empty files for msmc2
    kept_files = [f for f in os.listdir(out_dir) if os.path.isfile(os.path.join(out_dir, f)) and f.endswith("_msmc_input.txt") and os.path.getsize(os.path.join(out_dir, f)) > 0]

    if len(kept_files) == 0:
        print("Error! There are no usable file for msmc2, all inputs are empty!")
        exit(0)

    # Process each non-empty file
    for filename in kept_files:
        lines_to_keep = []
        with open(os.path.join(out_dir, filename), 'r') as file:
            for line in file:
                columns = line.strip().split('\t')
                if len(columns) >= 3 and int(columns[2]) > 0:
                    lines_to_keep.append(line)

        # Write back to the file with non-positive lines removed
        with open(os.path.join(out_dir, filename), 'w') as file:
            file.writelines(lines_to_keep)

        # Print information or run MSMC2 if needed
        if lines_to_keep:
            print(f"Processed {filename} - Non-positive lines removed.")
        else:
            print(f"{filename} has only non-positive lines and is now empty.")

    cmd5 = " ".join(["msmc2_Linux", ' '.join([os.path.join(out_dir, f) for f in kept_files]), "-o", out_dir+popid+"_msmc2"])

    print(cmd5)
    with open(popid+"_msmc2.log", 'w') as log:
        subprocess.run(cmd5, stdout=log, shell=True)

    # # Run MSMC2 combining information from all contigs
    # cmd5 = " ".join(["msmc2_Linux", out_dir+"*_msmc_input.txt -o", out_dir+popid+"_msmc2"])
    # print(cmd5)
    # with open(popid+"_msmc2.log", 'w') as log:
    #     p = subprocess.Popen(cmd5, stdout=log, shell=True)
    #     p.wait()

def psmc(ref_genome, contigs, popid, pop_ind, vcf, out_dir, mu, gen_time, p="10+22*2+4+6", num_cpus=4):
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
    "psmc_plot.pl", "-g", gen_time, "-x", "10", "-u", mu, "-M '"+",".join(pop_ind)+"'", popid, out_dir+'/'+popid+"_combined.psmc.final" ,";",
    "mv", popid+".*.par", out_dir, "; mv", popid+"*.eps", out_dir])
    print(cmd3)
    with open(out_dir+popid+"_psmc_combine.log", 'w') as log:
        p=subprocess.Popen(cmd3,stdout=log, shell=True)
    p.wait()

def smcpp(contigs, popid, pop_ind, vcf, out_dir, mu, gen_time):
    POP = popid+":"+",".join(pop_ind)
    if len(contigs) == 0:
        print("Error! No contigs to use! Make sure the threshold matches your data.")
        with open(out_dir+popid+"_model.final.json", 'w') as output:
            output.write("There was an error with SMC++. Please check the logs.")
        exit(0)

    for contig in contigs:
        cmd1 = ["smc++", "vcf2smc", vcf, out_dir+popid+"_"+contig+".smc.gz", contig, POP]
        print("LOG", out_dir+popid+"_"+contig+"_vcf2smc.log")
        with open(out_dir+popid+"_"+contig+"_vcf2smc.log", 'w') as log:
            p=subprocess.Popen(cmd1,stdout=log)
    p.wait()
    cmd2 = "".join(["smc++ estimate -o ", out_dir, popid, "_model.final.json ", mu, " ", out_dir, popid, "_*.smc.gz"])
    output2 = subprocess.check_output([cmd2], shell=True)
    output2 = output2.decode(encoding="utf-8")
    print(output2)
    cmd3 = "".join(["smc++ plot ", out_dir, popid, "_inference.png ", out_dir, popid, "_model.final.json/model.final.json -g", str(gen_time)," -c "])
    os.system(cmd3)
