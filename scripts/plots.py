#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot  the plots
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
import pandas as pd
import os
import re
from matplotlib.backends.backend_pdf import PdfPages

def plot_sfs(sfs, plot_title, output_file):
    plt.bar(np.arange(1, len(sfs)+1), sfs)
    plt.title(plot_title)
    plt.xlabel('# of alternate alleles')
    plt.ylabel('# sites')
    plt.xticks(np.arange(1, len(sfs)+1, 1.0))
    plt.savefig(output_file)
    plt.close()

def plot_stairwayplot2(popid, summary_file, out_dir):
    Ne_med=[]
    Ne1=[]
    Ne3=[]
    T=[]
    with open(summary_file) as input_file:
        line = input_file.readline()
        line = input_file.readline()
        while line != '':
            #if float(line.split('\t')[6]) not in Ne or float(line.split('\t')[5]) not in T:
            Ne_med.append(float(line.split('\t')[6]))
            Ne1.append(float(line.split('\t')[7]))
            Ne3.append(float(line.split('\t')[8]))
            T.append(float(line.split('\t')[5]))
            line = input_file.readline()
    plt.plot(T,Ne_med,color="blue")
    plt.plot(T,Ne1,color="grey")
    plt.plot(T,Ne3,color="grey")
    plt.title(popid)
    plt.xlabel('Time (in years)')
    plt.ylabel('Ne')
    plt.savefig(out_dir+popid+"_stw_plot.png")
    plt.close()

def plot_msmc2(popid, summary_file, mu, gen_time, out_dir):
    #scaled_left_bound = np.array()
    mu = float(mu)
    gen_time = float(gen_time)
    df = pd.read_csv(summary_file, delimiter='\t')
    # Scale Time
    # from units of the per-generation mutation rate to generations
    # to do this, divid by mutation rate
    scaled_left_bound = df['left_time_boundary']/mu * gen_time
    scaled_right_bound = df['right_time_boundary']/mu * gen_time
    # Pop Size
    scaled_pop = 1/df['lambda']
    pop_size_change = scaled_pop / (2*mu)
    plt.plot(scaled_left_bound, pop_size_change,color="blue")
    plt.title(popid)
    plt.xlabel('Time (in years)')
    plt.ylabel('Ne')
    plt.savefig(out_dir+popid+"_msmc2_plot.png")
    plt.close()

def plot_psmc(popid, sample_names, psmc_output_file,
              plot_output_prefix, gen_time, mut_rate,
              out_dir, x=0):
    # x: minimum generations, 0 for auto [10000]
    # R: DO not remove temporary files
    cmd = " ".join(["psmc_plot.pl -g", str(gen_time), "-x", str(x),
                    "-u", str(mut_rate), "-R", plot_output_prefix, psmc_output_file])
    os.system(cmd)
    os.system(f"cp {plot_output_prefix}.eps {out_dir+popid}_psmc_plot.eps")
    
def plot_distrib_gq(popid, gq, out_dir_gq):
    gq = OrderedDict(sorted(gq.items()))
    names = list(gq.keys())
    values = list(gq.values())
    plt.figure(figsize = [10, 5])
    plt.bar(range(len(gq)), values, tick_label=names)
    plt.title(popid + " GQ distribution")
    plt.xlabel('GQ')
    plt.ylabel('Numbre of sites')
    plt.savefig(out_dir_gq+popid + "_gq_distrib.png")
    plt.close()

def genotyping_coverage_plot(popid, snp_coverage, out_dir_stats, nb_plots=None, filter_prefix=None, nb_bins=None):
    """
    Generate a multi-page PDF report containing mean variants coverage plot for (selected) chromosomes.

    Parameters:
    - popid (str): Identifier for the population or sample.
    - snp_coverage (dict): Dictionary containing genotyping coverage data.
                          Format: {chromosome: {start_position: mean_coverage}}
    - out_dir_stats (str): Output directory path where the PDF report will be saved.
    - nb_plots (int, optional): Number of top chromosomes to include in the report. Default is None (include all).
    - filter_prefix (str, optional): Regular expression to filter chromosomes by name. Default is None (no filtering).
    - nb_bins (int, optional): Number of bins to group data points for better visualization. Default is None (no binning).

    Returns:
    None

    Saves a multi-page PDF report with coverage plots for the specified genotyping data.
    """
    if not os.path.exists(out_dir_stats):
        os.makedirs(out_dir_stats)

    # Sort contigs by CHR name
    sorted_contigs = sorted(snp_coverage.keys())

    # Filter contigs based on nb_plots and filter_prefix
    if nb_plots is not None:
        sorted_contigs = sorted_contigs[:nb_plots]

    if filter_prefix is not None:
        sorted_contigs = [chrm for chrm in sorted_contigs if re.match(filter_prefix, chrm)]

    chrm_count = len(sorted_contigs)
    plots_per_page = 8
    pages_needed = chrm_count // plots_per_page + int(chrm_count % plots_per_page > 0)

    with PdfPages(out_dir_stats + popid + "_variants_dist.pdf") as pdf:
        for page in range(pages_needed):
            fig, axes = plt.subplots(4, 2, figsize=(15, 5 * min(plots_per_page, chrm_count - page * plots_per_page)))

            for i in range(plots_per_page):
                chrm_idx = page * plots_per_page + i
                if chrm_idx >= chrm_count:
                    break

                chrm = sorted_contigs[chrm_idx]

                x_values = []
                y_values = []

                for start_pos_segment in snp_coverage[chrm]:
                    mean_coverage_for_this_segment = snp_coverage[chrm][start_pos_segment]
                    y = mean_coverage_for_this_segment
                    x = start_pos_segment
                    x_values.append(x)
                    y_values.append(y)

                # Calculate bin size based on nb_bins
                if nb_bins is not None:
                    bin_size = len(x_values) // nb_bins + 1
                    x_values = np.array(x_values)
                    y_values = np.array(y_values)

                    x_values_binned = [x_values[i:i + bin_size].mean() for i in range(0, len(x_values), bin_size)]
                    y_values_binned = [y_values[i:i + bin_size].mean() for i in range(0, len(y_values), bin_size)]

                    x_values = x_values_binned
                    y_values = y_values_binned
                
                axes[i // 2, i % 2].plot(x_values, y_values, label="Chromosome " + chrm)
                axes[i // 2, i % 2].set_xlabel("Start Position of Segment")
                # axes[i // 2, i % 2].set_ylabel("Mean Variants Prop")
                # axes[i // 2, i % 2].set_title("Chromosome " + chrm + " Average Variant Prop. per Segment")
                axes[i // 2, i % 2].set_ylabel("Dist. between SNPs")
                axes[i // 2, i % 2].set_title("Chromosome " + chrm + " distance between variants")

            # Adjust the height of the figure for a more standard paper ratio
            fig.set_size_inches(15, 8)

            plt.tight_layout()
            pdf.savefig()
            plt.close()
            
def plot_smcpp(popid, summary_file, out_dir):
    Ne=[]
    T=[]
    with open(summary_file) as input_file:
        line = input_file.readline()
        line = input_file.readline()
        while line != '':
            #if float(line.split('\t')[6]) not in Ne or float(line.split('\t')[5]) not in T:
            Ne.append(float(line.split(',')[2]))
            T.append(float(line.split(',')[1]))
            line = input_file.readline()
    plt.plot(T,Ne,color="blue")
    plt.title(popid)
    plt.xlabel('Time (in years)')
    plt.ylabel('Ne')
    plt.savefig(out_dir+popid+"_smcpp_plot.png")
    plt.close()
    
def plot_dadi_output_three_epochs(dadi_vals_list,name_pop,out_dir, mu, L, gen_time,
                                  xlim = None, ylim = None,
                                  max_v = -10**6, nb_plots_max = 10, title="Dadi pop. estimates"):
    best_model = None
    scenarios = {}
    for elem in dadi_vals_list:
        # elem is structured like this:
        # [iter_number, log_likelihood, [nuB, nuF, TB, TF], theta]
        log_lik = elem[1]
        # store scenarios sorted by their log_lik
        scenarios[float(log_lik)] = elem
    best_scenarios = sorted(scenarios.keys(), reverse = True)[:nb_plots_max]
    lines_x = []
    lines_y = []
    for scenario in best_scenarios:
        elem = scenarios[scenario]
        log_lik = elem[1]
        #popt = np.array(elem[2])
        theta = elem[3]
        Nanc = theta / (4*mu*L)
        # rescale to Nancestral effective size
        nuB = np.array(elem[2][0]) * Nanc
        nuF = np.array(elem[2][1]) * Nanc
        # Convert to years
        TB = np.array(elem[2][2]) * gen_time * Nanc
        TF = np.array(elem[2][3]) * gen_time * Nanc
        # now create x and y
        lines_x.append([0, TF, TF, TF+TB, TF+TB, (TF+TB)+0.01])
        lines_y.append([nuF, nuF, nuB, nuB, 1, 1])
    for i in range(1, len(lines_x)):
        x = lines_x[i]
        y = lines_y[i]
        plt.plot(x, y, '--', alpha = 0.4)
    #best model :
    best_x = lines_x[0]
    best_y = lines_y[0]
    plt.plot(best_x, best_y, 'r-', lw=2)
    plt.ylabel("Individuals (Na)")
    plt.xlabel("Time (years)")
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    plt.title(title)
    plt.savefig(out_dir+"/"+name_pop+"_dadi_plot.png")
    plt.close()

def Gplot(T_scaled_gen,gen_time,dadi_vals_list,name_pop,out_dir,popid, summary_file,summary_file2, max_v = -10**6, title="estimates",):
    import dadi
    with open(T_scaled_gen) as fo:
        line = fo.readline()
        fo.close()
    t_scaled_gen= line[1:-1]
    T_scaled_years=float(str.split(t_scaled_gen,", ")[-1])

    for dadi_params in dadi_vals_list:
        dadi_params[2][2]*=T_scaled_years
        dadi_params[2][3]*=T_scaled_years

    best_model = None
    for elem in dadi_vals_list:
        popt = elem[2]
        if elem[1] > max_v:
            max_v = elem[1]
            best_model = elem[2]
    best_x = [0, best_model[3], best_model[3], best_model[3]+best_model[2], best_model[3]+best_model[2], (best_model[3]+best_model[2])+0.01]
    best_y = [best_model[1], best_model[1], best_model[0], best_model[0], 1, 1]
    #with open(T_scaled_gen) as fo:
    #    line = fo.readline()
    #    fo.close()
    #t_scaled_gen= line[1:-1]
    #t_scaled_gen=str.split(t_scaled_gen,", ")[-1]
    Nanc=str.split(t_scaled_gen,", ")[0]
    #t_scaled_gen=[i*float(t_scaled_gen) for i in best_x]
    plt.plot(best_x, best_y,color="green")
    Ne=[]
    T=[]
    with open(summary_file2) as input_file2:
        line = input_file2.readline()
        line = input_file2.readline()
        while line != '':
            #if float(line.split('\t')[6]) not in Ne or float(line.split('\t')[5]) not in T:
            Ne.append(float(line.split(',')[2]))
            T.append(float(line.split(',')[1]))
            line = input_file2.readline()
    plt.plot(T,[i/(2*Ne[-1]) for i in Ne],color="blue")

    Ne_med=[]
    Ne1=[]
    Ne3=[]
    T=[]
    with open(summary_file) as input_file:
        line = input_file.readline()
        line = input_file.readline()
        while line != '':
            #if float(line.split('\t')[6]) not in Ne or float(line.split('\t')[5]) not in T:
            Ne_med.append(float(line.split('\t')[6]))
            Ne1.append(float(line.split('\t')[7]))
            Ne3.append(float(line.split('\t')[8]))
            T.append(float(line.split('\t')[5]))
            line = input_file.readline()
    plt.plot(T,[i/(2*Ne_med[-1]) for i in Ne_med],color="red")
    #plt.plot(T,[i/(2*Ne1[-1]) for i in Ne1],color="grey")
    #plt.plot(T,[i/(2*Ne3[-1]) for i in Ne3],color="grey")
    plt.title(popid)
    #plt.xlim(0,18000)
    plt.xlabel('Time (in years)')
    plt.ylabel('Na')
    plt.savefig(out_dir+popid+"_all_plot.png")
    plt.close()

def plot_pca(plink_eigenvec, plink_eigenval, popid, out_dir):
    # Load eigenvectors and eigenvalues
    with open(plink_eigenvec) as input:
        for line in input:
            if line.startswith("#"):
                continue
            else:
                num_cols = len(line.split())
                break
    eigenvectors = np.genfromtxt(plink_eigenvec, skip_header=1, usecols=range(2, num_cols))  # Skip header and first two columns
    eigenvalues = np.loadtxt(plink_eigenval)

    # Determine the number of components
    num_components = eigenvectors.shape[1]

    # Calculate the percentage of variance explained by each PC
    variance_explained = (eigenvalues / np.sum(eigenvalues)) * 100

    # Plot the first two principal components
    plt.figure(figsize=(8, 6))
    plt.scatter(eigenvectors[:, 0], eigenvectors[:, 1], c='blue', marker='o', s=50)
    plt.title(f'PCA: PC1 vs PC2 ({num_components} components)')
    plt.xlabel(f'PC1 ({variance_explained[0]:.2f}%)')
    plt.ylabel(f'PC2 ({variance_explained[1]:.2f}%)')
    plt.savefig(out_dir+"/"+popid+"_PCA.png")

    # Bar plot of explained variance
    plt.figure(figsize=(10, 6))
    plt.bar(range(1, num_components + 1), variance_explained * 100, color='blue', alpha=0.7)
    plt.xlabel('Number of Components (K)')
    plt.ylabel('Proportion of Explained Variance (%)')
    plt.title(f'Explained Variance by Components ({num_components} components)')
    plt.savefig(out_dir+"/"+popid+"_explained_var.png")
