"""
plots.py - Module of DemInfHelper for generating various plots and visualizations as part of the demography inference workflow.

This module provides functions for generating plots and visualizations related to population genetics, demographic inference,
and data analysis. It includes functions for Principal Component Analysis (PCA), Dadi model parameter estimation plots,
genotype quality distribution plots, and more.

Functions:
    - plot_pca(plink_eigenvec, plink_eigenval, popid, out_dir, n_clusters=9):
        Perform PCA on genetic data and generate plots, including PCA scatter plots and explained variance bar plots.
    
    - plot_dadi_output_three_epochs(dadi_vals_list, name_pop, out_dir, mu, L, gen_time,
                                   xlim=None, ylim=None, max_v=-10**6, nb_plots_max=10, title="Dadi pop. estimates"):
        Generate plots from Dadi model parameter estimates for three epochs, including parameter value distributions and trajectories.
    
    - plot_psmc(popid, sample_names, psmc_output_file, plot_output_prefix, gen_time, mut_rate, out_dir, x=0):
        Generate PSMC (Pairwise Sequentially Markovian Coalescent) plots based on PSMC output files.
    
    - plot_msmc2(popid, summary_file, mu, gen_time, out_dir):
        Generate MSMC (Multiple Sequentially Markovian Coalescent) plots based on MSMC2 summary files.
    
    - plot_stairwayplot2(popid, summary_file, out_dir):
        Generate stairway plot-style demographic history plots based on summary files.

    - plot_smcpp(popid, summary_file, out_dir):
        Generate plots of effective population size (Ne) over time based on smc++ summary files.
    
    - plot_distrib_gq(popid, gq, out_dir_gq):
        Generate Genotyping Quality (GQ) distribution plots based on GQ values from VCF files.

    - plot_sfs(sfs, plot_title, output_file):
        Generate Site Frequency Spectrum (SFS) plots based on SFS data.

    - plot_dadi_output_three_epochs(dadi_vals_list, name_pop, out_dir, mu, L, gen_time,
                                   xlim=None, ylim=None, max_v=-10**6, nb_plots_max=10, title="Dadi pop. estimates"):
        Generate plots from Dadi model parameter estimates for three epochs, including parameter value distributions and trajectories.

Each function is designed to provide visualization capabilities for different aspects of the demography inference workflow.
"""


import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
import pandas as pd
import os
import re
from matplotlib.backends.backend_pdf import PdfPages
import plotly.express as px
from sklearn.cluster import KMeans
from scipy.spatial import ConvexHull
import plotly.graph_objects as go

def plot_straight_x_y(x,y):
    x_1 = [x[0]]
    y_1 = []
    for i in range(0, len(y)-1):
        x_1.append(x[i])
        x_1.append(x[i])
        y_1.append(y[i])
        y_1.append(y[i])
    y_1 = y_1+[y[-1],y[-1]]
    x_1.append(x[-1])
    return x_1, y_1

def barplot_sfs(sfs, output_file, xlab = "Allele frequency", ylab = "SNP counts",
                folded=True, title = "SFS", transformed = False, normalized = False):
    sfs_val = []
    n = len(sfs)
    sum_sites = sum(sfs)
    for k, ksi in enumerate(sfs):
        # k does not starts at 0 but 1
        k = k+1
        if transformed:
            ylab = r'$ \phi_i $'
            if folded:
                val = ((k*(2*n - k)) / (2*n))*(ksi)
            else:
                val = ksi * k
        else:
            val = ksi
        sfs_val.append(val)

        if not transformed and not normalized:
            ylab = r'$ \eta_i $'
       
    #terminal case, same for folded or unfolded
    if transformed:
        last_bin = sfs[n-1] * n/2
    else:
        last_bin = sfs[n-1]
    sfs_val[-1] = last_bin
    if normalized:
        ylab = r'$ \phi_i $'
        sum_val = sum(sfs_val)
        for k, sfs_bin in enumerate(sfs_val):
            sfs_val[k] = sfs_bin / sum_val

    title = title+" (n="+str(len(sfs_val))+") [folded="+str(folded)+"]"+" [transformed="+str(transformed)+"]"
    print("SFS =", sfs)
    if folded:
        xlab = "Minor allele frequency"
    if transformed:
        print("Transformed SFS ( n =",len(sfs_val), ") :", sfs_val)
        #plt.axhline(y=1/n, color='r', linestyle='-')
    else:
        if normalized:
            # then plot a theoritical distribution as 1/i
            expected_y = [1/(2*x+1) for x in list(sfs.keys())]
            print(sum(expected_y))
            #plt.plot([x for x in list(sfs.keys())], expected_y, color='r', linestyle='-')
            #print(expected_y)
            
    x = [x+1 for x in range(len(sfs))]
    y = sfs_val
    # plt.xticks(x)
    plt.bar(x, y)
    plt.ylabel(ylab)
    plt.xlabel(xlab)
    plt.title(title)
    plt.savefig(output_file)
    plt.close()
    #plt.show()


def plot_sfs(sfs, plot_title, output_file):
    """
    Plot the Site Frequency Spectrum (SFS).

    This function generates a bar plot of the Site Frequency Spectrum (SFS) based on the input data.

    Parameters:
        sfs (list): A list of integers representing the SFS. Each value corresponds to the number of alternate alleles
                    observed at different sites.
        plot_title (str): The title for the plot.
        output_file (str): The path to save the generated plot as an image file.

    Example:
        sfs = [10, 15, 8, 3, 2]
        plot_title = "SFS Plot"
        output_file = "sfs_plot.png"
        plot_sfs(sfs, plot_title, output_file)

    Note:
        - The 'sfs' list should contain counts of alternate alleles at different sites.
        - The 'plot_title' will be used as the title of the generated plot.
        - The 'output_file' should specify the path to save the plot image.

    Returns:
        None: The function saves the plot as an image file but does not return any values.
    """
    plt.bar(np.arange(1, len(sfs)+1), sfs)
    plt.title(plot_title)
    plt.xlabel('# of alternate alleles')
    plt.ylabel('# sites')
    plt.xticks(np.arange(1, len(sfs)+1, 1.0))
    plt.savefig(output_file)
    plt.close()

def plot_stairwayplot2(popid, summary_file, out_dir, xlog=True, ylog=True, plot_format="png"):
    """
    Generate a stairway plot from summary data.

    This function reads summary data from a Stairway Plot analysis, extracts the
    effective population size (Ne) estimates over time, and generates a plot.

    Parameters:
        popid (str): Identifier for the population.
        summary_file (str): Path to the Stairway Plot summary data file.
        out_dir (str): Directory where the stairway plot image will be saved.

    Note:
        - The plot will be saved as "{popid}_stw_plot.png" in the specified 'out_dir'.

    Returns:
        None: The function saves the stairway plot as an image file but does not return any values.
    """
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
    plt.plot(T,Ne_med,color="red")
    plt.plot(T,Ne1,color="grey")
    plt.plot(T,Ne3,color="grey")
    plt.title(fr"[StairwayPlot2] {popid}")
    if xlog:
        plt.xscale("log")
    if ylog:
        plt.yscale("log")
    plt.xlabel('Time (years ago)')
    plt.ylabel('Ne')
    plt.savefig(out_dir+popid+f"_stw_plot.{plot_format}")
    plt.close()

def plot_msmc2(popid, summary_file, mu, gen_time, out_dir, xlog=True, ylog=True, plot_format="png"):
    """
    Generate a plot from MSMC2 summary data.

    This function reads summary data from MSMC2 analysis, scales the time and
    population size estimates, and generates a plot showing effective population size (Ne)
    changes over time.

    Parameters:
        popid (str): Identifier for the population.
        summary_file (str): Path to the MSMC2 summary data file.
        mu (float): Mutation rate per generation.
        gen_time (float): Generation time in years.
        out_dir (str): Directory where the plot image will be saved.

    Note:
        - The 'summary_file' should contain tab-separated summary data with columns
          'left_time_boundary', 'right_time_boundary', and 'lambda'.
        - Time is scaled from units of the per-generation mutation rate to generations.
        - The plot will be saved as "popid_msmc2_plot.png" in the specified 'out_dir'.

    Returns:
        None: The function saves the plot as an image file but does not return any values.
    """
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
    
    scaled_left_bound, pop_size_change = plot_straight_x_y(list(scaled_left_bound), list(pop_size_change))

    plt.plot(scaled_left_bound, pop_size_change,color="red")
    if xlog:
        plt.xscale("log")
    if ylog:
        plt.yscale("log")
    plt.title(fr"[MSMC2] {popid} $\mu$={mu} gen.t={gen_time}")
    plt.xlabel('Time (years ago)')
    plt.ylabel('Ne')
    plt.savefig(out_dir+popid+f"_msmc2_plot.{plot_format}")
    plt.close()

def plot_psmc(popid, sample_names, psmc_output_file,
              plot_output_prefix, gen_time, mut_rate,
              out_dir, kwargs):
    """
    Generate a PSMC plot from PSMC output.

    This function generates a PSMC plot using the PSMC plot script provided by PSMC
    (psmc_plot.pl). It takes PSMC output data and creates a plot in both EPS and PNG formats.

    Parameters:
        popid (str): Identifier for the population.
        sample_names (list): List of sample names used in PSMC analysis.
        psmc_output_file (str): Path to the PSMC output file.
        plot_output_prefix (str): Prefix for the plot output files.
        gen_time (float): Generation time in years.
        mut_rate (float): Mutation rate per generation.
        out_dir (str): Directory where the plot image files will be saved.
        x (int, optional): Minimum generations (0 for auto, default is 0).

    Note:
        - The 'sample_names' list should contain the names of samples used in PSMC analysis.
        - The 'psmc_output_file' is the output file generated by PSMC.
        - The 'plot_output_prefix' is used as the prefix for the output plot files.
        - The function assumes that 'psmc_plot.pl' is available in the system's PATH.

    Returns:
        None: The function generates PSMC plot files but does not return any values.
    """
    # x: minimum generations, 0 for auto [10000]
    # R: DO not remove temporary files
    cmd = " ".join(["psmc_plot.pl -g", str(gen_time), str(kwargs),
                    "-u", str(mut_rate), "-R", plot_output_prefix, psmc_output_file])
    os.system(cmd)
    print(cmd)
    os.system(f"cp {plot_output_prefix}.eps {out_dir+popid}_psmc_plot.eps")
    
def plot_distrib_gq(popid, gq, out_dir_gq, plot_format="png"):
    """
    Plot the distribution of Genotype Quality (GQ) values from a VCF file.

    Genotype Quality (GQ) is a phred score that represents the probability that the retained genotype was called
    correctly. It mirrors the efficiency of the calling process, which is influenced by external factors such as
    sequencing depth and coverage.

    This function takes a dictionary of GQ values and creates a bar plot showing the distribution of GQ values along
    with the corresponding count of sites. The GQ values provide insight into the quality of genotyping in the VCF file.

    Parameters:
        popid (str): Identifier for the population.
        gq (dict): Dictionary containing GQ values as keys and their counts as values.
        out_dir_gq (str): Directory where the GQ distribution plot will be saved.

    Note:
        - The 'gq' dictionary should contain GQ values as keys and the number of sites with each GQ value as values.
        - The function creates a bar plot showing the GQ distribution.
        - The plot is saved as "popid_gq_distrib.png" in the specified 'out_dir_gq'.

    Returns:
        None: The function generates the GQ distribution plot but does not return any values.
    """

    gq = OrderedDict(sorted(gq.items()))
    names = list(gq.keys())
    values = list(gq.values())
    plt.figure(figsize = [10, 5])
    plt.bar(range(len(gq)), values, tick_label=names)
    plt.title(f"{popid} GQ distribution")
    plt.xlabel("GQ")
    plt.ylabel('Numbre of sites')
    plt.savefig(f"{out_dir_gq}/{popid}_GQ_distrib.{plot_format}")
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
            
def plot_smcpp(popid, summary_file, out_dir, xlog=True, ylog=True, plot_format="png"):
    """
    Plot the effective population size (Ne) trajectory using SMC++ summary output.

    SMC++ is a tool for estimating demographic history from whole-genome sequencing data. This function reads
    SMC++ summary output files and creates a plot to visualize the effective population size (Ne) over time.

    Parameters:
        popid (str): Identifier for the population.
        summary_file (str): Path to the .csv SMC++ summary output file containing Ne and time information.
        out_dir (str): Directory where the SMC++ plot will be saved.

    Note:
        - The summary_file should be in the format generated by SMC++ containing Ne and time information.
        - The function creates a line plot showing the Ne trajectory over time.
        - The plot is saved as "popid_smcpp_plot.png" in the specified 'out_dir'.

    Returns:
        None: The function generates the SMC++ plot but does not return any values.
    """
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
    T, Ne = plot_straight_x_y(T,Ne)
    plt.plot(T,Ne,color="red")
    if xlog:
        plt.xscale("log")
    if ylog:
        plt.yscale("log")
    plt.title(fr"[SMC++] {popid}")
    plt.xlabel('Time (years ago)')
    plt.ylabel('Ne')
    plt.savefig(out_dir+popid+f"_smcpp_plot.{plot_format}")
    plt.close()
    
def plot_dadi_output_three_epochs(dadi_vals_list,name_pop,out_dir, mu, L, gen_time,
                                  xlim = None, ylim = None, xlog = True, ylog = True,
                                  max_v = -10**6, nb_plots_max = 10, title=None, plot_format="png"):
    """
    Plot demographic scenarios estimated using Dadi for a population with three epochs of size changes.

    Dadi is a tool for inferring demographic history from genetic data. This function takes a list of Dadi output values,
    extracts demographic parameters, and creates a plot to visualize population size changes over time.

    Parameters:
        dadi_vals_list (list): List of Dadi output values, where each element is structured as [iter_number, log_likelihood, [nuB, nuF, TB, TF], theta].
        name_pop (str): Name or identifier for the population.
        out_dir (str): Directory where the Dadi plot will be saved.
        mu (float): Mutation rate per site per generation.
        L (int): Total number of sites.
        gen_time (float): Generation time in years.
        xlim (tuple): Tuple specifying the x-axis limits (optional).
        ylim (tuple): Tuple specifying the y-axis limits (optional).
        max_v (int): Maximum value for the log-likelihood (optional, default is -10^6).
        nb_plots_max (int): Maximum number of scenarios to plot (optional, default is 10).
        title (str): Title for the Dadi plot (optional, default is "Dadi pop. estimates").

    Note:
        - The dadi_vals_list should contain Dadi output values, where each element is structured as described.
        - This function is specific to demographic scenarios with three epochs of population size changes.
        - The plot shows population size (individuals) over time (years) for different scenarios.
        - The best-scoring scenario is shown in red, and other scenarios are shown in dashed lines.
        - The plot is saved as "name_pop_dadi_plot.png" in the specified 'out_dir'.

    Returns:
        None: The function generates the Dadi plot but does not return any values.
    """
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
        lines_x.append([0, TF, TF, TF+TB, TF+TB])
        lines_y.append([nuF, nuF, nuB, nuB, 1])
    for i in range(1, len(lines_x)):
        x = lines_x[i]
        y = lines_y[i]
        plt.plot(x, y, '--', alpha = 0.4)
        # draws an artificial dotted line to represent infinity of ancestral pop
        plt.hlines(y[-1], x[-1], x[-1]*1.05, linestyle=':', alpha = 0.4)
    #best model :
    best_x = lines_x[0]
    best_y = lines_y[0]
    plt.plot(best_x, best_y, 'r-', lw=2)
    # draws an artificial dotted line to represent infinity of ancestral pop
    plt.hlines(best_y[-1], best_x[-1], best_x[-1]*1.05, linestyle=':', color='red')
    if xlog:
        plt.xscale("log")
    if ylog:
        plt.yscale("log")
    plt.ylabel("Individuals (Na)")
    plt.xlabel("Time (years ago)")
    if title:
        plt.title(title)
    else:
        plt.title(fr"[Dadi] {name_pop} $\mu$={mu} gen.t={gen_time}")
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    plt.savefig(out_dir+"/"+name_pop+f"_dadi_plot.{plot_format}")
    plt.close()

def Gplot(gen_time, mu, L, out_dir, popid, name_pop=None,
          dadi_file=None, summary_file=None, summary_file2=None,
          msmc2_summary_file=None, psmc_output_prefix=None,
          max_v=-10**6, title="estimates", plot_format="png"):
    """
    Generate a combined demographic inference plot overlaying all available methods.

    Each method is optional: the function checks whether the corresponding output file
    exists and silently skips any method whose file is absent. At least one method must
    produce data, otherwise the plot is not saved.

    Methods:
        dadi            -- requires dadi_file (.InferDM.bestfits)
        SMC++           -- requires summary_file2 (CSV)
        StairwayPlot2   -- requires summary_file (*final.summary)
        MSMC2           -- requires msmc2_summary_file (*_msmc2.final.txt)
        PSMC            -- requires psmc_output_prefix (*_combined.psmc.final)

    All plotted curves are extended to the oldest time point across methods with a
    dashed horizontal line. StairwayPlot2 and PSMC are also extended to t=0.

    Parameters:
        gen_time (float): Generation time in years.
        mu (float): Mutation rate per site per generation.
        L (float): Effective sequence length in bp (used to scale dadi Nanc).
        out_dir (str): Directory where the output PNG will be saved.
        popid (str): Population identifier used in the plot title and output filename.
        name_pop (str, optional): Unused, kept for call-site compatibility.
        dadi_file (str, optional): Path to the dadi bestfits file (.InferDM.bestfits).
        summary_file (str, optional): Path to the StairwayPlot2 final summary file
            (*final.summary, tab-separated, columns 5=T_years, 6=Ne_median).
        summary_file2 (str, optional): Path to the SMC++ CSV output file
            (columns 1=T_years, 2=Ne).
        msmc2_summary_file (str, optional): Path to the MSMC2 final txt file
            (tab-separated, columns left_time_boundary, lambda).
        psmc_output_prefix (str, optional): Path prefix for PSMC output; expects
            prefix + '_combined.psmc.final'. Parses the last iteration of each
            sample and averages across samples.
        max_v (float): Unused legacy parameter.
        title (str): Unused legacy parameter (plot title comes from popid).

    Returns:
        None: Saves the plot as out_dir + popid + '_combined_plot.png', or prints a warning
        and returns early if no method output files are found.
    """
    mu = float(mu)
    L = float(L)
    gen_time = float(gen_time)
    plt.figure(figsize=(10, 7))

    plotted = {}

    # --- dadi ---
    if dadi_file and os.path.exists(dadi_file):
        from parsing import dadi_output_parse
        dadi_vals_list = dadi_output_parse(dadi_file)
        if dadi_vals_list:
            scenarios = {float(e[1]): e for e in dadi_vals_list}
            best_elem = scenarios[sorted(scenarios.keys(), reverse=True)[0]]
            theta = best_elem[3]
            Nanc = theta / (4 * mu * L)
            TB = float(np.array(best_elem[2][2]) * gen_time * Nanc)
            TF = float(np.array(best_elem[2][3]) * gen_time * Nanc)
            T_dadi_max = TF + TB
            nuB = float(np.array(best_elem[2][0]) * Nanc)
            nuF = float(np.array(best_elem[2][1]) * Nanc)
            plt.plot([0, TF, TF, T_dadi_max, T_dadi_max], [nuF, nuF, nuB, nuB, Nanc],
                     color="green", label=r"$\partial$a$\partial$i")
            plotted['dadi'] = (T_dadi_max, Nanc)

    # --- SMC++ ---
    if summary_file2 and os.path.exists(summary_file2):
        Ne_smcpp = []
        T_smcpp = []
        with open(summary_file2) as f:
            f.readline()
            line = f.readline()
            while line != '':
                Ne_smcpp.append(float(line.split(',')[2]))
                T_smcpp.append(float(line.split(',')[1]))
                line = f.readline()
        if T_smcpp:
            T_smcpp, Ne_smcpp = plot_straight_x_y(T_smcpp, Ne_smcpp)
            plt.plot(T_smcpp, Ne_smcpp, color="blue", label="SMC++")
            plotted['smcpp'] = (T_smcpp[-1], Ne_smcpp[-1])

    # --- StairwayPlot2 ---
    if summary_file and os.path.exists(summary_file):
        Ne_med = []
        T_stw = []
        with open(summary_file) as f:
            f.readline()
            line = f.readline()
            while line != '':
                Ne_med.append(float(line.split('\t')[6]))
                T_stw.append(float(line.split('\t')[5]))
                line = f.readline()
        if T_stw:
            plt.plot(T_stw, Ne_med, color="red", label="StairwayPlot2")
            plt.hlines(Ne_med[0], 0, T_stw[0], colors="red", linestyles="dashed")
            plotted['stw'] = (T_stw[-1], Ne_med[-1])

    # --- MSMC2 ---
    if msmc2_summary_file and os.path.exists(msmc2_summary_file):
        df = pd.read_csv(msmc2_summary_file, delimiter='\t')
        T_msmc2 = list(df['left_time_boundary'] / mu * gen_time)
        Ne_msmc2 = list((1 / df['lambda']) / (2 * mu))
        if T_msmc2:
            T_msmc2, Ne_msmc2 = plot_straight_x_y(T_msmc2, Ne_msmc2)
            plt.plot(T_msmc2, Ne_msmc2, color="purple", label="MSMC2")
            plotted['msmc2'] = (T_msmc2[-1], Ne_msmc2[-1])

    # --- PSMC --- parse .psmc file directly, average across samples
    # Uses same scaling as psmc_plot.pl (opts_s = fq2psmcfa -s default = 100):
    #   N0 = theta / opts_s / (4 * mu)
    #   T  = 2 * N0 * t_left * gen_time
    #   Ne = lambda * N0
    if psmc_output_prefix:
        psmc_file = psmc_output_prefix + "_combined.psmc.final"
        if os.path.exists(psmc_file):
            opts_s = 100  # fq2psmcfa default window size in bp (-s flag)
            with open(psmc_file) as f:
                content = f.read()
            blocks = [b.strip() for b in content.split('//') if b.strip()]
            parsed_blocks = []
            for block in blocks:
                n_rd = None; theta = None; rs_data = []
                for line in block.strip().split('\n'):
                    parts = line.split('\t')
                    tag = parts[0]
                    if tag == 'RD':
                        n_rd = int(parts[1])
                    elif tag == 'TR':
                        theta = float(parts[1])
                    elif tag == 'RS':
                        rs_data.append((float(parts[2]), float(parts[3])))
                if n_rd is not None and theta is not None and rs_data:
                    parsed_blocks.append({'n_rd': n_rd, 'theta': theta, 'rs': rs_data})
            if parsed_blocks:
                max_rd = max(b['n_rd'] for b in parsed_blocks)
                best_blocks = [b for b in parsed_blocks if b['n_rd'] == max_rd]
                all_curves = []
                for b in best_blocks:
                    N0 = b['theta'] / opts_s / (4 * mu)
                    T_vals = [2 * N0 * t * gen_time for t, lam in b['rs']]
                    Ne_vals = [lam * N0 for t, lam in b['rs']]
                    pairs = [(t, n) for t, n in zip(T_vals, Ne_vals) if t > 0]
                    if pairs:
                        all_curves.append(([p[0] for p in pairs], [p[1] for p in pairs]))
                if all_curves:
                    if len(all_curves) == 1:
                        T_psmc_plot, Ne_psmc_plot = all_curves[0]
                    else:
                        T_common = np.linspace(
                            max(c[0][0] for c in all_curves),
                            min(c[0][-1] for c in all_curves),
                            500
                        )
                        interp_Ne = np.array([
                            np.interp(T_common, c[0], c[1]) for c in all_curves
                        ])
                        T_psmc_plot = list(T_common)
                        Ne_psmc_plot = list(interp_Ne.mean(axis=0))
                    T_psmc_plot, Ne_psmc_plot = plot_straight_x_y(T_psmc_plot, Ne_psmc_plot)
                    plt.plot(T_psmc_plot, Ne_psmc_plot, color="orange", label="PSMC (avg)")
                    plt.hlines(Ne_psmc_plot[0], 0, T_psmc_plot[0], colors="orange", linestyles="dashed")
                    plotted['psmc'] = (T_psmc_plot[-1], Ne_psmc_plot[-1])

    if not plotted:
        print(f"Warning: Combined plot for {popid} — no method output files found, skipping plot.")
        plt.close()
        return

    # --- dashed extensions to the most ancient time across all plotted methods ---
    T_max = max(v[0] for v in plotted.values())
    if 'dadi' in plotted:
        plt.hlines(plotted['dadi'][1], plotted['dadi'][0], T_max, colors="green", linestyles="dashed")
    if 'smcpp' in plotted:
        plt.hlines(plotted['smcpp'][1], plotted['smcpp'][0], T_max, colors="blue", linestyles="dashed")
    if 'stw' in plotted:
        plt.hlines(plotted['stw'][1], plotted['stw'][0], T_max, colors="red", linestyles="dashed")
    if 'msmc2' in plotted:
        plt.hlines(plotted['msmc2'][1], plotted['msmc2'][0], T_max, colors="purple", linestyles="dashed")
    if 'psmc' in plotted:
        plt.hlines(plotted['psmc'][1], plotted['psmc'][0], T_max, colors="orange", linestyles="dashed")

    plt.xscale("log")
    plt.yscale("log")
    plt.title(popid)
    plt.xlabel('Time (years ago)')
    plt.ylabel('Effective population size (Ne)')
    plt.legend()
    plt.savefig(out_dir + popid + f"_combined_plot.{plot_format}")
    plt.close()

def plot_pca(plink_eigenvec, plink_eigenval, popid, out_dir, n_clusters=9, plot_format="png"):
    """
    Perform Principal Component Analysis (PCA) on genetic data and generate plots and cluster information.

    This function performs PCA on genetic data using PLINK eigenvec and eigenval files, plots the first two principal components,
    performs k-means clustering, and generates informative plots and cluster information.

    Parameters:
        plink_eigenvec (str): Path to the PLINK eigenvec file containing PCA results.
        plink_eigenval (str): Path to the PLINK eigenval file containing eigenvalues.
        popid (str): Identifier for the population or dataset.
        out_dir (str): Directory where output files and plots will be saved.
        n_clusters (int): Number of clusters for k-means clustering (default is 9).

    Returns:
        None: The function generates various plots and cluster information but does not return any values.

    Note:
        - The function assumes that PLINK eigenvec and eigenval files are generated as part of PCA analysis.
        - It performs k-means clustering on the first two principal components.
        - Cluster information is saved in a CSV file named "popid_kn_clusters.csv".
        - Various plots, including PCA scatter plots and explained variance bar plot, are saved in the 'out_dir'.
    """
    if n_clusters is None:
        raise ValueError("Error: Please set a number of clusters with --n_clust_kmeans.")
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

    # Extract sample names from the eigenvec file
    with open(plink_eigenvec) as input:
        sample_names = [line.split()[1] for line in input if not line.startswith("#")]
    
    # Calculate the percentage of variance explained by each PC
    variance_explained = (eigenvalues / np.sum(eigenvalues)) * 100

    # Create a DataFrame for plotting
    pca_df = pd.DataFrame({
        'PC1': eigenvectors[:, 0],
        'PC2': eigenvectors[:, 1],
        'Sample': sample_names,
    })

    # Perform k-means clustering
    print(f"Performing kmeans clustering with k={n_clusters} clusters...")
    kmeans = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)
    pca_df['Cluster'] = kmeans.fit_predict(eigenvectors[:, :2])  # Using only the first two PCs for clustering

    # Group the data by cluster
    cluster_groups = pca_df.groupby('Cluster')

    # Iterate over each cluster
    print("Sorted by cluster:")
    for cluster_id, cluster_data in cluster_groups:
        print(f"Cluster {cluster_id}:")

        # Iterate over each sample in the cluster
        for sample, row in cluster_data.iterrows():
            print(f"Sample: {row['Sample']}, Cluster: {row['Cluster']}")

        print()  # Add an empty line to separate clusters

    # Iterate over each sample
    print("Sorted by sample:")
    # Open a CSV file for writing
    with open(f"{out_dir}/{popid}_k{n_clusters}_clusters.csv", 'w') as csv_file:
        csv_file.write('Sample,Cluster\n')  # Write the header
        # Iterate over each sample
        for sample, row in pca_df.iterrows():
            sample_name = row['Sample']
            cluster = row['Cluster']
            print(f"Sample: {sample_name}, Cluster: {cluster}")
            csv_file.write(f'{sample_name},{cluster}\n')  # Write each sample and its cluster to the file

    # Plot the first two principal components
    plt.figure(figsize=(8, 6))
    plt.scatter(eigenvectors[:, 0], eigenvectors[:, 1], c='blue', marker='o', s=50)
    plt.title(f'PCA: PC1 vs PC2 ({num_components} components)')
    plt.xlabel(f'PC1 ({variance_explained[0]:.2f}%)')
    plt.ylabel(f'PC2 ({variance_explained[1]:.2f}%)')
    plt.savefig(out_dir+"/"+popid+f"_PCA.{plot_format}")

    # Create a Plotly scatter plot with different symbols for each cluster
    symbol_sequence = range(1,10)  # Define symbol sequence

    # Map cluster number to symbol
    pca_df['Symbol'] = pca_df['Cluster'].apply(lambda x: symbol_sequence[x % len(symbol_sequence)])

    fig = px.scatter(pca_df, x='PC1', y='PC2', color='Cluster', symbol='Symbol',
                    title=f'PCA with k={n_clusters} Clusters: PC1 vs PC2 ({num_components} components)',
                    labels={'PC1': f'PC1 ({variance_explained[0]:.2f}%)', 'PC2': f'PC2 ({variance_explained[1]:.2f}%)'},
                    template='plotly_white', text='Sample')  # Add 'text' parameter for sample names

    # Add convex hulls
    for i in range(n_clusters):
        group = pca_df[pca_df['Cluster'] == i]
        points = group[['PC1', 'PC2']].values

        # Check if there are at least 3 points to form a convex hull
        if len(points) > 2:
            hull = ConvexHull(points)
            for simplex in hull.simplices:
                fig.add_trace(go.Scatter(x=points[simplex, 0], y=points[simplex, 1], mode='lines',
                                         line=dict(color=px.colors.qualitative.Set1[i % len(px.colors.qualitative.Set1)]),
                                         showlegend=False))
        else:
            print(f"Cluster {i} has fewer than 3 points; skipping convex hull.")

    # Update the layout
    fig.update_traces(marker=dict(size=10), selector=dict(type='scatter'))
    fig.update_layout(coloraxis_showscale=False)
    fig.update_layout(showlegend=False)
    # Save the interactive HTML plot
    fig.write_html(out_dir+"/"+popid+"_PCA_interactive.html")

    # Bar plot of explained variance
    plt.figure(figsize=(10, 6))
    plt.bar(range(1, num_components + 1), variance_explained * 100, color='blue', alpha=0.7)
    plt.xlabel('Number of Components (K)')
    plt.ylabel('Proportion of Explained Variance (%)')
    plt.title(f'Explained Variance by Components ({num_components} components)')
    plt.savefig(out_dir+"/"+popid+f"_explained_var.{plot_format}")
