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
from matplotlib.pyplot import stem
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

def build_defaults(out_dir, name_pop, program_path=None):
    """
    Build default derived paths from out_dir and name_pop.

    All output subdirectories and tool-specific paths follow a fixed layout rooted
    at out_dir. Users only need to set out_dir and name_pop in their config file;
    all other paths are inferred automatically unless explicitly overridden.

    Parameters:
        out_dir (str): Top-level output directory (e.g. './my_pop/').
        name_pop (list): List of population identifiers.
        program_path (str, optional): Path to the deminfhelper installation directory.
            Defaults to the directory containing this file.

    Returns:
        dict: Default values for all derived config keys.
    """
    if program_path is None:
        program_path = os.path.dirname(os.path.abspath(__file__))
    base = out_dir.rstrip('/') + '/'
    p0 = name_pop[0]  # used for single-pop default paths
    return {
        'out_dir_sfs':          base + 'output_sfs/',
        'path_to_sfs':          base + 'output_sfs/SFS_' + p0 + '.fs',
        'out_dir_stairwayplot2': base + 'output_stairwayplot2/',
        'summary_file_stw':     base + 'output_stairwayplot2/' + p0 + '/' + p0 + '.final.summary',
        'out_dir_dadi':         base + 'output_dadi/',
        'out_dir_msmc2':        base + 'output_msmc2/',
        'out_dir_smcpp':        base + 'output_smcpp/',
        'out_dir_psmc':         base + 'output_psmc/',
        'plot_file_smcpp':      base + 'output_smcpp/' + p0 + '_inference.csv',
        'out_dir_gq_distrib':   base + 'output_stats/',
        'final_out_dir':        base + 'inferences/',
        'out_dir_stats':        base + 'output_stats/',
        'path_to_stairwayplot2': program_path + '/bin/stairway_plot_es/',
        'blueprint_template':   program_path + '/bin/template.blueprint',
    }


def parse_config(config_file, args=None):
    """
    Parse a configuration file and return the configuration parameters as a dictionary.

    This function reads a configuration file line by line and extracts parameter names and values.
    The resulting dictionary contains configuration parameters with their corresponding values.

    Parameters:
    - config_file (str): The path to the configuration file.
    - args (dict): The dictionnary of command line arguments. 

    Returns:
    - param (dict): A dictionary containing configuration parameters and values.

    """

    param = {}
    with open(config_file, "rt") as config:
        for line in config:
            # if the line can be stripped, it's not empty
            if line.strip():
                if line[0] != "#":
                    param[line[:-1].split(": ")[0]] = line[:-1].split(": ")[1].strip()

    param["folded"]=bool(param["folded"])
    #param["transformed"]=bool(param["transformed"])
    param["name_pop"] = param["name_pop"].split(",")
    param["npop"]=int(param["npop"])
    if "out_dir" not in param:
        param["out_dir"] = "./" + param["name_pop"][0].strip() + "/"
    if "n_clust_kmeans" in param and param["n_clust_kmeans"] != None:
        param["n_clust_kmeans"] = eval(param["n_clust_kmeans"])
    if "missingness_by_sample" in param:
        param["missingness_by_sample"] = float(param["missingness_by_sample"])
    else:
        param["missingness_by_sample"] = 1.0
    if "missingness_by_site" in param:
        param["missingness_by_site"] = int(param["missingness_by_site"])
    else:
        param["missingness_by_site"] = 0
    if "cpus" in param:
        param["cpus"]=int(param["cpus"])
    else:
        param["cpus"]=None
    if "percentile_cutoff" not in param:
        param["percentile_cutoff"] = args.percentile_cutoff
    if "mask" not in param:
        if args:
            if args.mask != None:
                param["mask"]=args.mask
            else:
                param["mask"]=None
        else:
            param["mask"]=None
    param["sample_size"] = 0
    for p in param["name_pop"]:
        if p in list(param.keys()):
            param[p] = [item.strip() for item in param[p].split(",")]
        else:
            # if the pop is not defined in the config, the same list
            # of all samples from the VCF are used for every pop
            param[p] = get_sample_names(vcf=param["vcf"])
        param["n_"+p] = len(param[p])
        param["sample_size"] += param["n_" + p]


    # SETTING SOME DEFAULTS
    if "length_cutoff" not in param:
        # default contig size to keep is 1Mb
        param["length_cutoff"] = 100000
    else:
        param["length_cutoff"] = int(param["length_cutoff"])
    if "ref_genome" not in param:
        param["ref_genome"] = None
    if 'length_cutoff' not in param.keys():
        param["length_cutoff"] = 100000

    # Fill in any derived paths not explicitly set in the config
    for key, value in build_defaults(param['out_dir'], param['name_pop']).items():
        if key not in param:
            param[key] = value

    return param

def update_config(config_dict, config_file, args):
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

    initial_dict = {key.lower(): key for key in parse_config(config_file, args=args)}

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
            num_bins, mode, species_name = file.readline().strip().split()
            num_bins = int(num_bins)
            # Read the spectrum data
            spectrum_data = list(map(int, file.readline().strip().split()))
            # Check if the number of bins in the spectrum matches the expected number
            if len(spectrum_data) != num_bins:
                print("Len SFS=", len(spectrum_data), "Nb SFS bins. = ", num_bins)
                raise ValueError("Error: Number of bins in the spectrum doesn't match the specified number of bins (nb of haplotypes + 1).")
            # Read the mask data
            mask_data = list(map(int, file.readline().strip().split()))

            # Check if the size of the mask matches the number of bins in the spectrum
            if len(mask_data) != num_bins:
                raise ValueError("Error: Size of the mask doesn't match the number of bins in the spectrum.")
            # Apply the mask to the spectrum
            masked_spectrum = [spectrum_data[i] for i in range(num_bins) if not mask_data[i]]
    # Error handling
    except FileNotFoundError:
        print(f"Error: File not found - {sfs_file}")
    except ValueError as ve:
        print(f"Error: {ve}")
    except Exception as e:
        print(f"Error: {e}")
    # final return of SFS as a list
    return masked_spectrum

def get_contigs_lengths(vcf, length_cutoff=100000, only_names=False, contig_regex=None):
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
    with gzip.open(vcf, 'rt') as vcf_stream:
        line = vcf_stream.readline()
        while line != "":
            if line[0:8] == "##contig":
                contig_length = int(re.split('[=,]', line)[-1][:-2])
                contig_name = re.split('[=,]', line)[2]
                # keep only contigs that are longer than the length_cutoff parameter
                if contig_length >= length_cutoff:
                    contigs[contig_name] = contig_length
            elif line.startswith("#"):
                line = vcf_stream.readline()
                continue
            else:
                break
            line = vcf_stream.readline()
    #If not ##contig in the comments, need to estimate parsing all the file
    if len(contigs)==0:
        print("Warning: No ##contig info in VCF header. Need to parse the whole file...")
        pbar = tqdm(total=0, dynamic_ncols=True, unit='line', unit_scale=True) # Initialize the progress bar
        contigs_dict = {}
        with gzip.open(vcf, 'rt') as vcf_stream:
            line = vcf_stream.readline()
            while line != "":
                if line.startswith("#"):
                    line = vcf_stream.readline()
                    continue
                else:
                    line = line.split()
                    pos = line[1]
                    chrm = line[0]
                    contigs[chrm] = pos
                # at the end, change line
                pbar.update(1)
                line = vcf_stream.readline()
        pbar.close()
    print("Finished getting contig sizes.")
    if contig_regex:
        contig_regex = re.compile(contig_regex)
        print(f"Keeping only contigs with regex: {contig_regex}")
        contigs = {contig_name: contig_size for contig_name, contig_size in contigs.items() if (re.match(contig_regex, contig_name))}
    
    
    # print(f"Contigs={contigs}")
    if only_names:
        return list(contigs.keys())
    return contigs

def get_sample_names(vcf):
    """Get samples list from the VCF header
    """
    print(f"Parsing {vcf} to get samples list.")
    # Parsing VCF in gzip format
    with gzip.open(vcf, 'rt') as vcf:
        line = vcf.readline()
        while line != "":
            if line[0:6] == "#CHROM":
                samples = list(line.split())[9:]
                break
            line = vcf.readline()
    return samples

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
                 keep_modified_vcf = False, modified_vcf_ram = False, mem=4096):
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
                    " --memory ", str(mem), \
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


def parse_bed(bed_file):
    mask = {}
    with open(bed_file, "r") as bed_stream:
        for line in bed_stream:
            contig, start, end = line.split()[0:3]
            if contig not in mask.keys():
                mask[contig] = [int(start), int(end)]
            else:
                mask[contig] += [int(start), int(end)]
    return mask


def pos_in_mask(kept_pos, target_pos):
    if kept_pos is None:
        return True
    if target_pos not in kept_pos:
        return False
    return True


def kept_pos(mask, chrm):
    """Given a chromosome, will return a list of all the kept positions in the mask
    """
    kept_pos = []
    if chrm not in mask.keys():
        return kept_pos
    for k, pos in enumerate(mask[chrm][:-1:2]):
        kept_pos+=list(range(pos, mask[chrm][k+1]))
    return kept_pos


def missingness(genotypes, stats):
    nb_missing = 0
    for index, elem in enumerate(genotypes):
        if elem[1] == ":":
            g1, g2 = True, True  # if the genotype is missing, both alleles are considered missing
        else:
            g1, g2 = elem[0] == ".", elem[2] == "."
        nb_missing += (g1 + g2)
        # print(elem[0], elem[2], nb_missing)
        stats["missingness_by_sample"][index] += (g1 or g2)
    # print(genotypes, nb_missing)
    stats["missingness_by_site"][nb_missing] += 1
    return nb_missing


def is_indel_or_pluriallelic(reformed_line, stats = None): 
    """
    Check if the variant is an indel or pluriallelic
    This function checks the variant represented by reformed_line to determine if it is an indel
    (insertion or deletion) or a pluriallelic site (multiple alleles).

    Parameters:
    - reformed_line (list): A list representing a VCF line with variant information.

    Returns:
    - bool: True if the variant is an indel or pluriallelic, False otherwise.
    """
    if "," in reformed_line[4]:
        if stats is not None:
            stats["pluriall"] += 1
        return True  # pluriallelic site
    elif len(reformed_line[4]) > 1 or len(reformed_line[3]) > 1:
        if stats is not None:
            stats["pluriall"] += 1
        return True  # indel
    return False  # neither indel nor pluriallelic


def exess_het(field):
    pattern = re.compile(r"ExcessHet=([-+]?\d*\.\d+|\d+)")
    match = pattern.search(field)
    return float(match.group(1))


def first_parsing(PARAM, contigs, mask = None):
    cols_in_vcf = {}    # Dictionary to store column indices for each population
    exclude_pos = {}  # Dictionary to store positions to exclude for each chromosome
    snp_dist_list = []  # List to store SNP distances
    tab_size = 4 # size of the tab to keep the last 4 positions
    stats = {"above_snp_distrib_cutoff": 0,
            # pos masked by bed file
            "masked_pos": 0,
            "missingness_by_sample": [0] *  PARAM["sample_size"],
            "missingness_by_site": [0] * (PARAM["sample_size"] * 2 + 1),
            "missing_data": 0,
            "pluriall": 0,
            "indels": 0
            }
    nb_line_not_skip = 0
    current_chrm = None  # Variable to keep track of the current chromosome
    pbar = tqdm(total=0, dynamic_ncols=True, unit='line', unit_scale=True) # Initialize the progress bar
    with gzip.open(PARAM["vcf"],  mode='rt') as vcf:
        line = vcf.readline()   
        # we read all the lines of the vcf
        while line != "":
            if line[0:6] == "#CHROM":
                # parsing the header line
                header_list = line[:-1].split("\t")
                for p in PARAM["name_pop"]:
                    pos_ind_pop = []
                    for ind in PARAM[p]:
                        pos_ind_pop.append(header_list.index(ind))
                    cols_in_vcf["pos_"+p] = pos_ind_pop
                # skip the rest of the code for this line
                line = vcf.readline()
                continue
            if line.startswith("#"):
                 # ignore comments
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
                # if mask:
                #     pos_kept = kept_pos(mask, chrm)
                # else:
                #     pos_kept = None
                # initialize exclude_pos for this chrm
                exclude_pos[chrm] = set()
            split_line = line.split("\t")
            if is_indel_or_pluriallelic(split_line): # if snp is noeither snp nor biallelic then it is skip
                line = vcf.readline()
                continue
            reformed_line = split_line[0:9]
            for sample_pos_list in cols_in_vcf.values():
                for sample_pos in sample_pos_list:
                    reformed_line.append(split_line[sample_pos])
            skip_line = False
            if line[0] != "#":
                # comment
                skip_line = False     
            if not skip_line:
                nb_line_not_skip += 1
                nb_missing = missingness(reformed_line[9:], stats) # calculate missingness
                if nb_missing > PARAM["missingness_by_site"]: # site missingness threshold
                    stats["missing_data"] += 1
                    skip_line = True # skip line if too much missing data
            #à quoi ca sert
            # if not pos_in_mask(pos_kept, pos):
            #     stats["masked_pos"] += 1
            #     skip_line = True
            if not skip_line:
                snp_nb += 1
                split_line = line.split("\t")
                for p in PARAM["name_pop"]:
                    # do something with this snp
                    tab[snp_nb % tab_size] = pos
                    if snp_nb > tab_size:
                        snp_dist = tab[snp_nb % tab_size] - tab[(snp_nb+1) % tab_size]
                        snp_dist_list.append(snp_dist)
            else:
                # if the line is skipped, put it in the exclude list
                exclude_pos[chrm].add(pos)
            line = vcf.readline()
            pbar.update(1)
    pbar.close()  # Close the progress bar when done
    print("First parsing the VCF to determine the distribution of distances between variants.")
    stats['missingness_by_sample'] = np.array(stats['missingness_by_sample']) / nb_line_not_skip
    return cols_in_vcf, exclude_pos, snp_dist_list, stats, tab_size


import matplotlib.pyplot as plt

# Fonction 1 : Plot simple de "missing_sample"
def plot_missing_sample(data, PARAM):
    plt.figure(figsize=(6, 4))
    plt.hist(data, bins='auto', color='tomato', edgecolor='black')
    plt.title("Histogram of Missing Data Proportions (Sample)")
    plt.xlabel("Missing Proportion")
    plt.ylabel("Frequency")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    os.makedirs(PARAM ["out_dir_stats"], exist_ok=True)
    plt.savefig(PARAM ["out_dir_stats"] + "missing_sample_histogram.png")

# Fonction 2 : Plot cumulée en fréquence de "missing_site"
def plot_missing_site_cumulative(data, PARAM):
    cumulative = np.cumsum(data)
    cumulative_freq = cumulative / cumulative[-1]  # Normaliser
    plt.figure(figsize=(8, 4))
    plt.plot(range(len(data)), cumulative_freq, marker='o', linestyle='-', color='slateblue')
    plt.title("Missing Site - Cumulative Frequency (Original Order)")
    plt.xlabel("Index")
    plt.ylabel("Fréquence Cumulée")
    plt.grid(True)
    plt.tight_layout()
    os.makedirs(PARAM ["out_dir_stats"], exist_ok=True)
    plt.savefig(PARAM ["out_dir_stats"] + "missing_site_cumulative.png")

def delete_individuals(cols_in_vcf, PARAM, missingness_by_sample, missingness_threshold = 1.0):
    """
    Delete individuals from the VCF based on missingness threshold.

    This function checks the missingness of each individual in the VCF and removes those
    whose missingness exceeds the specified threshold.

    Parameters:
    - cols_in_vcf (dict): Dictionary containing column indices for each population.
    - PARAM (dict): Dictionary containing configuration parameters.
    - stats (dict): Dictionary containing statistics, including missingness by sample.
    - missingness_threshold (float): Threshold for missingness to determine if an individual should be removed.

    Returns:
    - None
    """
    # Calculate the total number of samples
    if missingness_by_sample is None:
        print("No missingness data provided. Skipping individual deletion.")
        return
    # Identify individuals to delete based on missingness threshold
    individuals_to_delete = []
    for i, missing in enumerate(missingness_by_sample):
        if missing  > missingness_threshold:
            individuals_to_delete.append(9 + i)
    sample_size = 0
    # Remove identified individuals from cols_in_vcf
    for pop in PARAM["name_pop"]:
        cols_in_vcf["pos_" + pop] = [i for i in cols_in_vcf["pos_" + pop] if i not in individuals_to_delete]
        PARAM["n_" + pop] = len(cols_in_vcf["pos_" + pop])
        sample_size += PARAM["n_" + pop]
    PARAM["sample_size"] = sample_size
    print(f"Deleted {len(individuals_to_delete)} individuals with missingness rate inferior to {missingness_threshold}.")


# Function using dynamic distance evaluation with rolling positions
def vcf_line_parsing(PARAM, SFS = False, GQ = False, SMCPP = False, mask = None, percentile_cutoff = 90):
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
    SFS_dict = {}
    GQ_dict = {}
    All_snp_count = 0
    Kept_snp_count = 0
    if mask:
        print(f"Omitting variants not present in {mask}.")
        mask = parse_bed(mask)
    # we initialize a sfs for each population
    for p in PARAM["name_pop"]:
        SFS_dict[p] = np.array(sfs.build_sfs(n=PARAM["sfs_size"], folded=PARAM["folded"], sfs_ini=True))
    if GQ:
        for p in PARAM["name_pop"]:
            GQ_dict[p] = {}
    cols_in_vcf, exclude_pos, snp_dist_list, stats, tab_size = first_parsing(PARAM, contigs, mask=mask)
    # delete_individuals(cols_in_vcf, PARAM, stats["missingness_by_sample"], PARAM["missingness_by_sample"])
   #def
    keeping_threshold = np.percentile(snp_dist_list, percentile_cutoff)
    print(f"SFS parsing: Done. Filtering variants keeping only variants with a distance that is lower than : {keeping_threshold} ({percentile_cutoff}' percentile)")
    snps_distance_by_chr = {}
    current_chrm = None
    pbar = tqdm(total=0, dynamic_ncols=True, unit='line', unit_scale=True) # Initialize the progress bar    
    with gzip.open(PARAM["vcf"],  mode='rt') as vcf:
        line = vcf.readline()    
        # we read all the lines of the vcf
        while line != "":
            split_line = line.split("\t")
            if line.startswith("#") or is_indel_or_pluriallelic(split_line): 
                # ignore other comments
                line = vcf.readline()
                continue
            # the line is a variant, count it
            # if float(split_line[7].split(";")[5].split("=")[1]) > 60.:
            #     line = vcf.readline()
            #     continue
            All_snp_count += 1
            chrm = line.split()[0]
            # if chrm not in contigs:
            #     line = vcf.readline()
            #     continue
            pos = int(split_line[1])
            if chrm != current_chrm:
                # new chromosome/contig initialization
                current_chrm = chrm
                snps_distance_by_chr[chrm] = {}
                # reset the tab, otherwise distances are going to be false
                tab = [1] * tab_size
                snp_nb = 0
            if pos not in exclude_pos[chrm]:
                snp_nb += 1
                for p in PARAM["name_pop"]:
                    # do something with this snp
                    tab[snp_nb % tab_size] = pos
                    if snp_nb > tab_size:
                        snp_dist = tab[snp_nb % tab_size] - tab[(snp_nb+1) % tab_size]
                        if snp_dist <= keeping_threshold:
                            Kept_snp_count += 1
                            snps_distance_by_chr[chrm][pos] = snp_dist
                            SFS_dict[p] = sfs.build_sfs(PARAM["sfs_size"], folded=PARAM["folded"],  sfs_ini = False, \
                                                        line = split_line, sfs = SFS_dict[p], pos_ind = cols_in_vcf["pos_"+p])
                        else:
                            stats["above_snp_distrib_cutoff"] += 1
                if GQ:
                    for p in PARAM["name_pop"]:
                        GQ_dict[p] = distrib_GQ(GQ_pop = GQ_dict[p], line = split_line, pos_ind = cols_in_vcf["pos_"+p])
            line = vcf.readline()
            pbar.update(1)
    pbar.close()  # Close the progress bar when done
    # print(SFS_dict)
    L = (All_snp_count - Kept_snp_count) / All_snp_count * Nseq
    print("Finished building SFS.")
    print(f"Stats:\n-------\n# SNPs:\t{All_snp_count}\nMissing sample:\t{stats['missingness_by_sample']}\nMissing site:\t{stats['missingness_by_site']}\nPluriallelic sites or indels:\t{stats['pluriall']}"+\
          f"\nSites above percentile cutoff [{percentile_cutoff}]:\t{stats['above_snp_distrib_cutoff']}\n"+\
          f"Kept sites:\t{Kept_snp_count}")
    plot_missing_sample(stats['missingness_by_sample'], PARAM)
    plot_missing_site_cumulative(stats['missingness_by_site'], PARAM)
    return SFS_dict, GQ_dict, round(L), snps_distance_by_chr
