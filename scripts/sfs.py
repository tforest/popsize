"""
Site Frequency Spectrum (SFS) Module of DemInfHelper.

This module provides functions for computing and transforming the Site Frequency Spectrum (SFS) from genetic data.
The SFS is a fundamental concept in population genetics, representing the distribution of allele frequencies within a population.

Functions:
- build_sfs(n, folded, sfs_ini, line=[], sfs=[], pos_ind=None): Compute the SFS based on genetic data.
- transform_sfs(sfs, n, folded): Transform the SFS, folded or unfolded.

Usage:
1. Import this module: `import sfs`
2. Use the functions to compute and manipulate SFS data.
"""
import itertools
import numpy as np

def build_sfs(n, folded, sfs_ini, line=[], sfs=[], pos_ind = None):
    """
    Compute the Site Frequency Spectrum (SFS).

    This function computes the Site Frequency Spectrum (SFS) based on the provided parameters and genetic data.
    It can either initialize an empty SFS or update an existing one based on the input VCF line and other parameters.

    Parameters:
    - n (int): The number of samples.
    - folded (bool): True if the SFS should be folded, False for unfolded.
    - sfs_ini (bool): True to initialize an empty SFS, False to update an existing one.
    - line (list): A list representing a VCF line split into fields.
    - sfs (list): A list containing the initial SFS.
    - pos_ind (list): Positions for each sample in the VCF.

    Returns:
    - sfs (list): The updated Site Frequency Spectrum.
    """
    if sfs_ini:
        if folded:
            return [0]+[0] * n
        else:
            return [0]+[0] * (2 * n - 1)
    else:
        geno_ind = [line[i] for i in pos_ind]
        gen = list(itertools.chain.from_iterable([[int(i[0]), int(i[2])] for i in geno_ind])) #we get the genotypes
        #gen a list of 1 and 0
        if folded:
            count=min(sum(gen), 2*n-sum(gen)) #if folded, we count the minor allele
        else:
            count=sum(gen) #if not folded, we count the alternate allele (1)
        # warning! SFS[0] are monomorphic sites, as well as SFS[len(gen)] is SFS is unfoled
        # if count != 0 and (folded) or (not folded and count != 2*n) # these are monomorphic sites
        sfs[count] += 1 #the sfs is incremented
        return(sfs)


def transform_sfs(sfs, n, folded):
    """
    Transform the Site Frequency Spectrum (SFS).

    This function transforms the Site Frequency Spectrum (SFS) using a specific formula, which depends on
    whether the SFS is folded or unfolded.

    Parameters:
    - sfs (list): The original Site Frequency Spectrum.
    - n (int): The number of samples.
    - folded (bool): True if the SFS is folded, False if unfolded.

    Returns:
    - sfs_trans (list): The transformed SFS.

    Note: The function applies a formula to transform the SFS and returns the transformed SFS.
    """
    #to transform the vcf
    if folded:
        sfs_trans = []
        for i in range(1,len(sfs)+1):
            if i<len(sfs):
                sfs_trans.append(round(i * (2 * n - i) * sfs[i-1] / (2 * n)))
            else:
                sfs_trans.append(round(i*sfs[i-1]))
        return sfs_trans
    else:
        return [round(i * sfs[i-1]) for i in range(len(sfs))]


