# -*- coding: utf-8 -*-
"""
parsing() : parses the vcd, compute the sfs if SFS=True, output the msmc input file
if MSMC=False
"""

import gzip
import os
import numpy as np
# from inferences import *

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


def distrib_GQ(GQ_pop, line = [], pos_ind = None): #PL is a dict
    samples = [line[i] for i in pos_ind]
    gq_line = [i[4:] for i in samples] #we get the genotypes
    gq_line = [i.split(",") for i in gq_line]
    for sublist in gq_line:
        gq2 = [int(i) for i in sublist]
        gq2 = [i - np.min(gq2) for i in gq2]
        gq2.remove(0)
        min = np.min(gq2)
        bin = min - int(str(min)[-1])
        if bin in GQ_pop.keys():
            GQ_pop[bin] = GQ_pop[bin]+1
        else:
            GQ_pop[bin] = 1
    return(GQ_pop)


def parse_config(config_file):
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
    Update the configuration file with the values from the given config_dict.
    Preserve in-place comment lines and add new entries at the end.

    Parameters:
    - config_dict (dict): The dictionary containing the updated configuration values.
    - config_file (str): The path to the configuration file.

    Returns:
    - None
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
    Parses a site frequency spectrum (SFS) file and returns a masked spectrum.

    Parameters:
    - sfs_file (str): The path to the SFS file to be parsed.

    Returns:
    - list: A masked spectrum, where bins indicated by the mask are excluded.

    Raises:
    - FileNotFoundError: If the specified SFS file is not found.
    - ValueError: If there are inconsistencies in the file format or data, such as mismatched bin counts.

    Example:
    >>> parse_sfs('example.fs')
    [1, 0, 3, 2, 0]
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
    """Parsing gzipped VCF header to get contig lengths from the ##contig comments
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
# Function using segments    
# def vcf_line_parsing(PARAM, SFS = False, GQ = False, SMCPP = False, segments_size = 1000):
#     # cutoff is the minimum size of each contig to be used
#     # required for SMC++, as it works for contigs > 100kb or 1mb
#     length_cutoff = int(PARAM["length_cutoff"])

#     genome_segments = {}
#     # the genome is going to be chopped in segments of a maximum size
#     # defined by segments_size parameter.
#     # Each segment will store its SFS.
#     # as well as the mean coverage of SNPs
#     # At the end, only segments with sufficient coverage will be kept

#      ### L: Estimated number of sequence genotyped.
#      # from GADMA
#      # Assume total length of sequence that was used for SFS building is equal to Nseq.
#      # From this data total number of X SNP’s were received.
#      # But not all of them were used for SFS: some of them (e.g. total number of Y SNP’s) were filtered out.
#      # Then we should count filtered SNP’s and take L value the following way:
#      # L = (X - Y) / X * Nseq
    
#     with gzip.open(PARAM["vcf"],  mode='rt') as vcf:
#         SFS_dict = {}
#         GQ_dict = {}
#         cols_in_vcf = {}
        
#         segments_mean_snps_freq = []
#         Nseq = 0
#         All_snp_count = 0
#         Kept_snp_count = 0    
#         # we initialize a sfs for each population
#         for p in PARAM["name_pop"]:
#             genome_segments[p] = {}            
#             SFS_dict[p] = np.array(sfs.build_sfs(n=PARAM["n_"+str(p)], folded=PARAM["folded"], sfs_ini=True))
#         if GQ:
#             for p in PARAM["name_pop"]:
#                         GQ_dict[p] = {}
#         if SMCPP:
#             contigs = []
#         line = vcf.readline()
#         pbar = tqdm(total=0, dynamic_ncols=True, unit='line', unit_scale=True) # Initialize the progress bar
#         # we read all the lines of the vcf
#         while line != "":
#             if line[0:6] == "#CHROM":
#                 # parsing the header line
#                 for p in PARAM["name_pop"]:
#                     pos_ind_pop = []
#                     for ind in PARAM[p]:
#                         pos_ind_pop.append(line[:-1].split("\t").index(ind))
#                     cols_in_vcf["pos_"+p] = pos_ind_pop
#                 # skip the rest of the code for this line
#                 line = vcf.readline()
#                 continue
#             if line.startswith("#"):
#                 # ignore other comments
#                 line = vcf.readline()
#                 continue
#             # the line is a variant, count it
#             All_snp_count += 1
#             chrm = line.split()[0]
#             pos = int(line.split()[1])
#             if chrm not in genome_segments[p].keys():
#                 genome_segments[p][chrm] = {}
#                 start_pos = pos
#                 end_pos = start_pos+segments_size
#                 # initiate SFS for the first segment of this CHR
#                 # each segment is defined by its startpos, it belongs to a chromosome, itself from a certain (sub)population
#                 genome_segments[p][chrm][start_pos] = {}
#                 genome_segments[p][chrm][start_pos]['sfs'] = sfs.build_sfs(n=PARAM["n_"+str(p)], folded=PARAM["folded"], sfs_ini=True)
#                 genome_segments[p][chrm][start_pos]['count_snps'] = 0
#                 genome_segments[p][chrm][start_pos]['snps_mean'] = 0
#             elif pos > end_pos:
#                 # when we switch segment, same CHR,
#                 # We first store previous segment values.

#                 # Proportion of SNPs of all sites scanned for this segment
#                 # when switching segment: pos-start_pos should be
#                 # equal to segment_size, except for the last segment which can be shorter.
#                 genome_segments[p][chrm][start_pos]['snps_mean'] = genome_segments[p][chrm][start_pos]['count_snps'] / (pos-start_pos)
#                 segments_mean_snps_freq.append(genome_segments[p][chrm][start_pos]['snps_mean'])
#                 # change segment
#                 start_pos = pos
#                 end_pos = start_pos+segments_size
#                 # initiate SFS for this new segment
#                 genome_segments[p][chrm][start_pos] = {}
#                 genome_segments[p][chrm][start_pos]['sfs'] = sfs.build_sfs(n=PARAM["n_"+str(p)], folded=PARAM["folded"], sfs_ini=True)
#                 genome_segments[p][chrm][start_pos]['count_snps'] = 0
#                 genome_segments[p][chrm][start_pos]['snps_mean'] = 0

#             if line[0] != "#" and ".:" not in line and "/." not in line and "./" not in line and ".|" not in line and "|." not in line and "," not in line.split("\t")[4]:    #we only keep the bi-allelique sites
#                 Kept_snp_count += 1
#                 split_line = line.split("\t")
#                 for p in PARAM["name_pop"]:
#                     genome_segments[p][chrm][start_pos]['sfs'] = sfs.build_sfs(n=PARAM["n_"+p], \
#                                                                                 folded=PARAM["folded"], \
#                                                                                 sfs_ini = False, \
#                                                                                 line = split_line, \
#                                                                                 sfs = genome_segments[p][chrm][start_pos]['sfs'], \
#                                                                                 pos_ind = cols_in_vcf["pos_"+p])
#                     genome_segments[p][chrm][start_pos]['count_snps'] += 1
#                     # SFS_dict[p] = sfs.build_sfs(n=PARAM["n_"+p], folded=PARAM["folded"],  sfs_ini = False, \
#                     #         line = split_line, sfs = SFS_dict[p], pos_ind = cols_in_vcf["pos_"+p])

#                 if GQ:
#                     for p in PARAM["name_pop"]:
#                         GQ_dict[p] = distrib_GQ(GQ_pop = GQ_dict[p], line = split_line, pos_ind = cols_in_vcf["pos_"+p])
#             line = vcf.readline()
#             pbar.update(1)
#     pbar.close()  # Close the progress bar when done
#     print("SFS parsing: Done. Filtering variants keeping only segments with sufficient coverage.")
#     sorted_data = np.sort(segments_mean_snps_freq)
#     data_median = np.median(sorted_data)
#     Kept_snp_count = 0
#     snps_coverage_by_chr = {}
#     for p in PARAM["name_pop"]:
#         for chrm in genome_segments[p]:
#             snps_coverage_by_chr[chrm] = {}
#             for start_pos_segment in genome_segments[p][chrm]:
#                 snps_coverage_by_chr[chrm][start_pos_segment] = 0
#                 if genome_segments[p][chrm][start_pos_segment]['snps_mean'] >= data_median:
#                     # keep this segment, add its SFS to the final SFS for this p
#                     Kept_snp_count += 1
#                     SFS_dict[p] += np.array(genome_segments[p][chrm][start_pos_segment]['sfs'])
#                     snps_coverage_by_chr[chrm][start_pos_segment] = genome_segments[p][chrm][start_pos_segment]['snps_mean']
#     L = (All_snp_count - Kept_snp_count) / All_snp_count * Nseq
    
#     return SFS_dict, GQ_dict, round(L), snps_coverage_by_chr
# Function using dynamic distance evaluation with rolling positions
def vcf_line_parsing(PARAM, SFS = False, GQ = False, SMCPP = False, segments_size = 1000):
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
