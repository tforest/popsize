# -*- coding: utf-8 -*-
"""
parsing.py - Module of DemInfHelper for parsing and processing VCF files and other data formats in the demography inference workflow.

This module provides functions for parsing and processing VCF files, calculating Site Frequency Spectra (SFS), Genotyping Quality (GQ) distributions,
and performing various data preprocessing tasks. It also includes functions for dynamic distance evaluation, SNP filtering, and more.

Functions:
    - distrib_GQ(GQ_pop, line=[], pos_ind=None, bin_size=10):
        Parse Genotyping Quality (GQ) information from a VCF file and calculate GQ distributions.
    
    - parse_config(config_file):
        Parse a configuration file and return the configuration parameters as a dictionary.
    
    - update_config(config_dict, config_file):
        Update a configuration file with values from a dictionary and preserve in-place comments.
    
    - parse_sfs(sfs_file):
        Parse a Site Frequency Spectrum (SFS) file and return a masked spectrum.
    
    - get_contigs_lengths(vcf, length_cutoff=100000, only_names=False):
        Parse contig lengths from a VCF file header and return a dictionary of contig names and lengths.
    
    - dadi_output_parse(dadi_output_file):
        Parse Dadi model output files and return the results as a list of values.
    
    - pca_from_vcf(popid, vcf_file, nb_samples, out_dir, ploidy=2, keep_modified_vcf=False, modified_vcf_ram=False):
        Perform Principal Component Analysis (PCA) on genetic data from a VCF file and generate PCA plots.
    
    - vcf_line_parsing(PARAM, SFS=False, GQ=False, SMCPP=False, segments_size=1000):
        Parse VCF lines, calculate SNP distances, and generate Site Frequency Spectra (SFS) and GQ distributions.

Each function in this module is designed to assist in parsing and processing genetic data and configuration files for the demography inference workflow.
"""

import gzip
import os
import numpy as np

if __package__ is None or __package__ == '':
    import inferences
    import sfs
    import plots
else:
    from . import inferences
    from . import sfs
    from . import plots

import re
from tqdm import tqdm  # Import tqdm for the progress bar


def distrib_GQ(GQ_pop, line = [], pos_ind = None, bin_size = 10): #PL is a dict
    """
    Calculate the distribution of Genotyping Quality (GQ) values from a VCF line.

    This function parses GQ values from the INFO field of a VCF line, groups them into bins,
    and counts the number of GQ values in each bin.

    Parameters:
    - GQ_pop (dict): A dictionary to store the GQ distribution.
    - line (list): A list representing a VCF line split into fields.
    - pos_ind (list): List of positions for each sample in the VCF.
    - bin_size (int): The size of bins for grouping GQ values.

    Returns:
    - GQ_pop (dict): Updated GQ distribution dictionary.

    Raises:
    - ValueError: If the GQ field is not found in the VCF FORMAT field.
   """

    format_field = line[8].split(":")  # Split the FORMAT field
    
    if "GQ" not in format_field:
        raise ValueError("GQ field not found in FORMAT field")

    gq_index = format_field.index("GQ")  # Find the position of GQ field in FORMAT

    samples = [line[i] for i in pos_ind]  # Extract sample-specific information

    gq_values = [sample.split(":")[gq_index] for sample in samples]  # Extract GQ values

    for gq_value in gq_values:
        gq_value = int(gq_value)  # Convert GQ value to integer
        # Group GQ values into bins of specifed size (e.g. bin_size=10: 0-9, 10-19, etc.)
        bin_value = gq_value - (gq_value % bin_size) 

        if bin_value in GQ_pop:
            GQ_pop[bin_value] += 1
        else:
            GQ_pop[bin_value] = 1

    return GQ_pop

def parse_config(config_file):
    """
    Parse a configuration file and return the configuration parameters as a dictionary.

    This function reads a configuration file line by line and extracts parameter names and values.
    The resulting dictionary contains configuration parameters with their corresponding values.

    Parameters:
    - config_file (str): The path to the configuration file.

    Returns:
    - param (dict): A dictionary containing configuration parameters and values.

    """

    param = {}
    with open(config_file, "rt") as config:
        line=config.readline()
        while line != "":
            if line[0] != "#":
                    param[line[:-1].split(": ")[0]] = line[:-1].split(": ")[1].strip()
            line = config.readline()

    param["folded"]=bool(param["folded"])
    #param["transformed"]=bool(param["transformed"])
    param["name_pop"] = param["name_pop"].split(",")
    param["npop"]=int(param["npop"])
    if "n_clust_kmeans" in param and param["n_clust_kmeans"] != None:
        param["n_clust_kmeans"] = eval(param["n_clust_kmeans"])
    if "cpus" in param:
        param["cpus"]=int(param["cpus"])
    else:
        param["cpus"]=None
    for p in param["name_pop"]:
        param[p] = [item.strip() for item in param[p].split(",")]
        param["n_"+p] = len(param[p])


    # SETTING SOME DEFAULTS
    if "length_cutoff" not in param:
        # default contig size to keep is 1Mb
        param["length_cutoff"] = 100000
    else:
        param["length_cutoff"] = int(param["length_cutoff"])
    if "ref_genome" not in param:
        param["ref_genome"] = None
    if 'length_cutoff' not in param.keys():
        param["length_cutoff"] = length_cutoff
    return param

def update_config(config_dict, config_file):
    """
    Update a configuration file with values from a dictionary and preserve in-place comments.

    This function takes a dictionary of configuration parameters and their values and updates
    an existing configuration file. It preserves any comments in the file and adds new entries
    at the end if necessary.

    Parameters:
    - config_dict (dict): A dictionary containing updated configuration values.
    - config_file (str): The path to the configuration file to be updated.

    Returns:
    - None

    Note: The function modifies the configuration file in place.
    """
    # Read the existing config file
    with open(config_file, 'r') as file:
        lines = file.readlines()

    initial_dict = {key.lower(): key for key in parse_config(config_file)} 

    # Create a mapping of lowercase keys to original keys
    key_mapping = {key.lower(): key for key in config_dict}

    # Update values in the lines based on config_dict or add new entries
    for i in range(len(lines)):
        line = lines[i].strip()
        if line.startswith('#') or not line:
            # Preserve in-place comment lines
            continue

        key, value = map(str.strip, line.split(':', 1))

        # Use the original key if it's found in the config_dict
        key_original = key_mapping.get(key.lower(), key)

        if key_original in config_dict:
            # If the key is found in the config_dict, update the value
            updated_value = config_dict[key_original]
            # Check if the value is a list and format it accordingly
            if isinstance(updated_value, list):
                updated_value_str = ', '.join(map(str, updated_value))
            else:
                updated_value_str = updated_value
            lines[i] = f"{key_original}: {updated_value_str}\n"

    # Add new entries from config_dict at the end of the file if they are not already present
    for key, value in config_dict.items():
        # Check if the key is not already present in the file
        if key.lower() not in initial_dict.keys():
            lines.append(f"{key}: {value}\n")

    # Write the updated lines back to the config file
    with open(config_file, 'w') as output:
        for line in lines:
            output.write(line)


def parse_sfs(sfs_file):
    """
    Parse a Site Frequency Spectrum (SFS) file and return a masked spectrum.

    This function reads an SFS file, extracts the spectrum data, and applies a mask to it.
    The mask excludes specific bins from the spectrum, resulting in a masked SFS.

    Parameters:
    - sfs_file (str): The path to the SFS file to be parsed, in dadi's .fs format.

    Returns:
    - masked_spectrum (list): A masked SFS as a list of integers.

    Raises:
    - FileNotFoundError: If the specified SFS file is not found.
    - ValueError: If there are inconsistencies in the file format or data.

    Note: The actual structure of the SFS file is based on dadi's fs format.
    """
    try:
        with open(sfs_file, 'r') as file:
            # Read the first line which contains information about the file
            num_individuals, mode, species_name = file.readline().strip().split()
            num_individuals = int(num_individuals)
            # Read the spectrum data
            spectrum_data = list(map(int, file.readline().strip().split()))
            # Check if the number of bins in the spectrum matches the expected number
            if len(spectrum_data) != num_individuals:
                raise ValueError("Error: Number of bins in the spectrum doesn't match the expected number of individuals.")
            # Read the mask data
            mask_data = list(map(int, file.readline().strip().split()))

            # Check if the size of the mask matches the number of bins in the spectrum
            if len(mask_data) != num_individuals:
                raise ValueError("Error: Size of the mask doesn't match the number of bins in the spectrum.")
            # Apply the mask to the spectrum
            masked_spectrum = [spectrum_data[i] for i in range(num_individuals) if not mask_data[i]]
    # Error handling
    except FileNotFoundError:
        print(f"Error: File not found - {sfs_file}")
    except ValueError as ve:
        print(f"Error: {ve}")
    except Exception as e:
        print(f"Error: {e}")
    # final return of SFS as a list
    return masked_spectrum

def get_contigs_lengths(vcf, length_cutoff=100000, only_names = False):
    """
    Parse contig lengths from a VCF file header and return a dictionary of contig names and lengths.

    This function reads the header of a VCF file and extracts contig information, including contig names
    and their corresponding lengths. It returns a dictionary with contig names as keys and lengths as values.

    Parameters:
    - vcf (str): The path to the VCF file.
    - length_cutoff (int): The minimum contig length to be included.
    - only_names (bool): If True, return a list of contig names only.

    Returns:
    - contigs (dict or list): A dictionary of contig names and lengths or a list of contig names.

    Note: The function can return a dictionary of contig names and lengths or a list of contig names.
    """
    contigs = {}
    print(f"Parsing {vcf} to get contigs sizes.")
    # Parsing VCF in gzip format
    with gzip.open(vcf, 'rt') as vcf:
        line = vcf.readline()
        while line != "":
            if line[0:8] == "##contig":
                contig_length = int(re.split('[=,]', line)[-1][:-2])
                contig_name = re.split('[=,]', line)[2]
                # keep only contigs that are longer than the length_cutoff parameter
                if contig_length >= length_cutoff:
                    contigs[contig_name] = contig_length
            elif line.startswith("#"):
                line = vcf.readline()
                continue
            else:
                break
            line = vcf.readline()
    #If not ##contig in the comments, need to estimate parsing all the file
    if len(contigs)==0:
        print("Warning: No ##contig info in VCF header. Need to parse the whole file...")
        pbar = tqdm(total=0, dynamic_ncols=True, unit='line', unit_scale=True) # Initialize the progress bar
        contigs_dict = {}
        with gzip.open(vcf, 'rt') as vcf:
            line = vcf.readline()
            while line != "":
                if line.startswith("#"):
                    continue
                else:
                    line = line.split()
                    pos = line[1]
                    chrm = line[0]
                    conitgs[chrm] = pos
                # at the end, change line
                pbar.update(1)
                line = vcf.readline()
        pbar.close()
    if only_names:
        return list(contigs.keys())
    print("Finished getting contig sizes.")
    return contigs

def dadi_output_parse(dadi_output_file):
    """
    Parse Dadi model output files and return the results as a list of values.

    This function parses output files generated by Dadi demographic modeling and extracts the results,
    including log likelihoods, parameter values, and other relevant information.

    Parameters:
    - dadi_output_file (str): The path to the Dadi output file.

    Returns:
    - all_vals (list): A list of parsed values containing log likelihoods, parameter values, etc.

    Note: The structure of Dadi output files can vary, and parsing may need to be adapted accordingly.
    """
    all_vals = []
    ite = 0
    converged = False
    with open(dadi_output_file) as dadi_output:
        for line in dadi_output:
            # Log(likelihood)       nuB     nuF     TB      TF      misid   theta
            if line.startswith("#"):
                # check if converged
                if "Converged" in line:
                    converged = True
                if converged and len(all_vals) > 0:
                    # all converged values have been parsed
                    # skip the top 100 results that are printed as the second series of values
                    return all_vals
                # then skip if comments
                continue
            ite += 1
            line_parsed = line.strip().split()
            logL = float(line_parsed[0])
            nuB = float(line_parsed[1])
            nuF = float(line_parsed[2])
            TB = float(line_parsed[3])
            TF = float(line_parsed[4])
            theta = float(line_parsed[5])
            all_vals.append([ite, logL, [nuB, nuF, TB, TF], theta])
    return all_vals

def pca_from_vcf(popid, vcf_file, nb_samples, out_dir, ploidy = 2,
                 keep_modified_vcf = False, modified_vcf_ram = False):
    """
    Perform Principal Component Analysis (PCA) on genetic data from a VCF file and generate PCA plots.

    This function conducts PCA on genetic data from a VCF file, generates PCA plots, and saves the results
    in the specified output directory.

    Parameters:
    - popid (str): Identifier for the population being analyzed.
    - vcf_file (str): The path to the VCF file containing genetic data.
    - nb_samples (int): The number of samples in the VCF file.
    - out_dir (str): The directory where PCA results and plots will be saved.
    - ploidy (int): The ploidy level for the genetic data (default is 2).
    - keep_modified_vcf (bool): If True, keep the modified VCF file; otherwise, it will be deleted.
    - modified_vcf_ram (bool): If True, use RAM for the modified VCF file (only relevant for large VCF files).

    Returns:
    - None

    """
    plink_out_dir = out_dir+"/plink/"
    if not os.path.exists(plink_out_dir):
        os.makedirs(plink_out_dir)
    # need to use bcftools to add IDs to replace the "." with unique IDs for each variant 
    cmd1 = "".join(["bcftools annotate --set-id +'%CHROM:%POS' ", \
                    vcf_file, " -Oz -o ", \
                    plink_out_dir+popid+"_IDs.vcf.gz"])
    cmd2 = "".join(["plink2 --vcf ", plink_out_dir+popid+"_IDs.vcf.gz", \
                    " --make-bed --allow-extra-chr --max-alleles ", str(ploidy), \
                    " --snps-only --out ", plink_out_dir+popid, " --freq"])
    print(cmd1)
    os.system(cmd1)
    print(cmd2)
    os.system(cmd2)
    if keep_modified_vcf == False:
        # remove temporary file
        os.remove(plink_out_dir+popid+"_IDs.vcf.gz")
    cmd3 = "".join(["plink2 --bfile ", plink_out_dir+popid, \
                    " --pca ", str(nb_samples-1), \
                    " --out ", plink_out_dir+popid+".pca --allow-extra-chr --read-freq ", \
                     plink_out_dir+popid+".afreq"])
    print(cmd3)
    os.system(cmd3)
    # Generate plot
    plots.plot_pca(plink_out_dir+popid+".pca.eigenvec", plink_out_dir+popid+".pca.eigenval", popid = popid, out_dir = out_dir)

# Function using dynamic distance evaluation with rolling positions
def vcf_line_parsing(PARAM, SFS = False, GQ = False, SMCPP = False):
    """
    Parse VCF lines, calculate SNP distances, and generate Site Frequency Spectra (SFS) and GQ distributions.

    This function parses VCF lines, calculates SNP distances between samples, and can generate Site Frequency Spectra (SFS)
    and Genotyping Quality (GQ) distributions depending on the specified parameters.

    Parameters:
    - PARAM (dict): A dictionary containing configuration parameters.
    - SFS (bool): If True, generate Site Frequency Spectra (SFS).
    - GQ (bool): If True, generate Genotyping Quality (GQ) distributions.
    - SMCPP (bool): If True, perform specific tasks for SMC++ input.

    Returns:
    - None

    Example:
    >>> params = {'folded': True, 'name_pop': ['pop1', 'pop2'], ...}
    >>> SFS = True
    >>> GQ = True
    >>> vcf_line_parsing(params, SFS, GQ)

    Note: The function performs various tasks based on the specified parameters.
    """
    # cutoff is the minimum size of each contig to be used
    # required for SMC++, as it works for contigs > 100kb or 1mb
    length_cutoff = int(PARAM["length_cutoff"])
    contigs = get_contigs_lengths(vcf=PARAM["vcf"], length_cutoff=length_cutoff)
    # total size of all contigs used
    Nseq = sum(list(contigs.values()))
    tab_size = 4
    percentile = 90
    SFS_dict = {}
    GQ_dict = {}
    cols_in_vcf = {}
    current_chrm = None
    All_snp_count = 0
    Kept_snp_count = 0
    snp_dist_list = []
    # we initialize a sfs for each population
    for p in PARAM["name_pop"]:
        SFS_dict[p] = np.array(sfs.build_sfs(n=PARAM["n_"+str(p)], folded=PARAM["folded"], sfs_ini=True))
    if GQ:
        for p in PARAM["name_pop"]:
            GQ_dict[p] = {}
    pbar = tqdm(total=0, dynamic_ncols=True, unit='line', unit_scale=True) # Initialize the progress bar
    # First parse
    print("First parsing the VCF to determine the distribution of distances between variants.")
    with gzip.open(PARAM["vcf"],  mode='rt') as vcf:
        line = vcf.readline()    
        # we read all the lines of the vcf
        while line != "":
            if line[0:6] == "#CHROM":
                # parsing the header line
                for p in PARAM["name_pop"]:
                    pos_ind_pop = []
                    for ind in PARAM[p]:
                        pos_ind_pop.append(line[:-1].split("\t").index(ind))
                    cols_in_vcf["pos_"+p] = pos_ind_pop
                # skip the rest of the code for this line
                line = vcf.readline()
                continue
            if line.startswith("#"):
                # ignore other comments
                line = vcf.readline()
                continue
            chrm = line.split()[0]
            pos = int(line.split()[1])
            if chrm != current_chrm:
                # new chromosome/contig initialization
                current_chrm = chrm
                # reset the tab, otherwise distances are going to be false
                tab = [1] * tab_size
                snp_nb = 0
            if line[0] != "#" and ".:" not in line and "/." not in line and "./" not in line and ".|" not in line and "|." not in line and "," not in line.split("\t")[4]:    #we only keep the bi-allelique sites
                snp_nb += 1
                split_line = line.split("\t")
                for p in PARAM["name_pop"]:
                    # do something with this snp
                    tab[snp_nb % tab_size] = pos
                    if snp_nb > tab_size:
                        snp_dist = tab[snp_nb % tab_size] - tab[(snp_nb+1) % tab_size]
                        snp_dist_list.append(snp_dist)
            line = vcf.readline()
            pbar.update(1)
    pbar.close()  # Close the progress bar when done
    keeping_threshold = np.percentile(snp_dist_list, percentile)
    print(f"SFS parsing: Done. Filtering variants keeping only variants with a distance that is lower than : {keeping_threshold} ({percentile}' percentile)")
    snps_distance_by_chr = {}
    pbar = tqdm(total=0, dynamic_ncols=True, unit='line', unit_scale=True) # Initialize the progress bar    
    with gzip.open(PARAM["vcf"],  mode='rt') as vcf:
        line = vcf.readline()    
        # we read all the lines of the vcf
        while line != "":
            if line[0:6] == "#CHROM":
                # parsing the header line
                for p in PARAM["name_pop"]:
                    pos_ind_pop = []
                    for ind in PARAM[p]:
                        pos_ind_pop.append(line[:-1].split("\t").index(ind))
                    cols_in_vcf["pos_"+p] = pos_ind_pop
                # skip the rest of the code for this line
                line = vcf.readline()
                continue
            if line.startswith("#"):
                # ignore other comments
                line = vcf.readline()
                continue
            # the line is a variant, count it
            All_snp_count += 1
            chrm = line.split()[0]
            pos = int(line.split()[1])
            if chrm != current_chrm:
                # new chromosome/contig initialization
                current_chrm = chrm
                snps_distance_by_chr[chrm] = {}
                # reset the tab, otherwise distances are going to be false
                tab = [1] * tab_size
                snp_nb = 0
            if line[0] != "#" and ".:" not in line and "/." not in line and "./" not in line and ".|" not in line and "|." not in line and "," not in line.split("\t")[4]:    #we only keep the bi-allelique sites
                snp_nb += 1
                split_line = line.split("\t")
                for p in PARAM["name_pop"]:
                    # do something with this snp
                    tab[snp_nb % tab_size] = pos
                    if snp_nb > tab_size:
                        snp_dist = tab[snp_nb % tab_size] - tab[(snp_nb+1) % tab_size]
                        if snp_dist < keeping_threshold:
                            Kept_snp_count += 1
                            snps_distance_by_chr[chrm][pos] = snp_dist
                            SFS_dict[p] = sfs.build_sfs(n=PARAM["n_"+p], folded=PARAM["folded"],  sfs_ini = False, \
                                                        line = split_line, sfs = SFS_dict[p], pos_ind = cols_in_vcf["pos_"+p])
                if GQ:
                    for p in PARAM["name_pop"]:
                        GQ_dict[p] = distrib_GQ(GQ_pop = GQ_dict[p], line = split_line, pos_ind = cols_in_vcf["pos_"+p])
            line = vcf.readline()
            pbar.update(1)
    pbar.close()  # Close the progress bar when done
    L = (All_snp_count - Kept_snp_count) / All_snp_count * Nseq
    print("Finished building SFS.")
    return SFS_dict, GQ_dict, round(L), snps_distance_by_chr
