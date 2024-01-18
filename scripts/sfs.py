"""
build_sfs (n, folded,sfs_ini, line='', sfs=[]) : compute the sfs
transform_sfs(sfs, n, folded) : transforms the sfs
"""

import itertools
import numpy as np

def build_sfs(n, folded, sfs_ini, line=[], sfs=[], pos_ind = None):
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
    #to transform the vcf
    #formula from Elise thesis
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


