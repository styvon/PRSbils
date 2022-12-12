#!/usr/bin/env python3

"""
usage: PRSbils.py [-h] --ref_dir REF_DIR --bim_prefix BIM_PREFIX --sst_file
                  SST_FILE --map_file MAP_FILE [--map_gcol MAP_GCOL] --n_gwas
                  N_GWAS --out_dir OUT_DIR [--a A] [--b B] [--c C] [--d D]
                  [--e E] [--f F] [--n_iter N_ITER]
                  [--chrom CHROM [CHROM ...]] [--beta_std [BETA_STD]]
                  [--thres THRES] [--ignore_blk [IGNORE_BLK]]
                  [--verbose [VERBOSE]] [--fixtau [FIXTAU]] [--flip [FLIP]]

PRSbils: Polygenic Risk Score with Bilevel Continuous Shrinkage

optional arguments:
  -h, --help            show this help message and exit
  --ref_dir REF_DIR     Full path (including folder name) to the directory
                        (ldblk_1kg_eur, ldblk_1kg_eas or ldblk_1kg_afr) that
                        contains information on the LD reference panel
                        (snpinfo_1kg_hm3 and ldblk_1kg_chr*.hdf5)
  --bim_prefix BIM_PREFIX
                        Full path and the prefix of the bim file
  --sst_file SST_FILE   Full path and the file name of the GWAS summary
                        statistics Summary statistics file must be tab
                        delimited with the following format (including the
                        header line): SNP A1 A2 BETA P Or: SNP A1 A2 OR P
  --map_file MAP_FILE   Full path and the file name of the mapping file from
                        SNPs to sets. Map file must have the MSigDB file
                        format : [CHROM] [SNPID] [SETNAME] [GROUP (optional)]
  --map_gcol MAP_GCOL   (Not yet implemented) Column id for group info in map_file. Default is None
                        (no group information)
  --n_gwas N_GWAS       Sample size of the GWAS
  --out_dir OUT_DIR     Output directory and output filename prefix of the
                        posterior effect size estimates
  --a A                 Parameter a for sigma. Default is 0.001
  --b B                 Parameter b in the gamma prior for sigma. Default is
                        0.001
  --c C                 Parameter c for group implementation. Default is 0.001
  --d D                 Parameter d for group implementation. Default is 0.001
  --e E                 Parameter e for lambda. Default is 1
  --f F                 Parameter f for lambda. Default is 0.5
  --n_iter N_ITER       Number of maximum iterations. Default is 1000
  --chrom CHROM [CHROM ...]
                        The chromosome on which the model is fitted. Default
                        is iterating through 22 autosomes
  --beta_std [BETA_STD]
                        Boolean. If True, return standardized posterior SNP
                        effect sizes. Default is False
  --thres THRES         Threshold value for convergence. Default is 0.001
  --ignore_blk [IGNORE_BLK]
                        Boolean. If false, elements categorized as same set
                        but in different ld will be treated as different sets.
                        Default is False
  --verbose [VERBOSE]   Boolean. Default is False
  --fixtau [FIXTAU]     Boolean. If false, will estimate tau in algorithm,
                        otherwise will fix at 1. Default is False
  --flip [FLIP]         Boolean. If false, signs of beta in summary statistics
                        will not be flipped when there's a switch in A1 and A2
                        alleles. Default is False
"""

import sys
import argparse
import time
import parseg
import vbclass_mcmc
import gigrnd

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        print('Invalid Boolean argument\n')
        sys.exit(2)

def parse_param():
    parser = argparse.ArgumentParser(description='PRSbils: Polygenic Risk Score with Bilevel Continuous Shrinkage')
    parser.add_argument("--ref_dir", help="Full path (including folder name) to the directory (ldblk_1kg_eur, ldblk_1kg_eas or ldblk_1kg_afr) that contains information on the LD reference panel (snpinfo_1kg_hm3 and ldblk_1kg_chr*.hdf5)", required=True)
    parser.add_argument("--bim_prefix", help="Full path and the prefix of the bim file", required=True)
    parser.add_argument("--sst_file",
        help='''
        Full path and the file name of the GWAS summary statistics
        Summary statistics file must be tab delimited with the following format (including the header line):

        SNP          A1   A2   BETA      P

        Or:
        SNP          A1   A2   OR        P
        ''', required=True)
    parser.add_argument("--map_file", help='''
        Full path and the file name of the mapping file from SNPs to sets.
        Map file must have the MSigDB file format :

        [CHROM]    [SNPID]    [SETNAME]    [GROUP (optional)]
        ''', required=True)
    parser.add_argument("--map_gcol", type=int, default=None, help="Column id for group info in map_file. Default is None (no group information)", required=False)
    parser.add_argument("--n_gwas", type=int, help="Sample size of the GWAS", required=True)
    parser.add_argument("--out_dir", help="Output directory and output filename prefix of the posterior effect size estimates", required=True)
    parser.add_argument("--a", type=float, default=0.001, help="Parameter a for sigma. Default is 0.001", required=False)
    parser.add_argument("--b", type=float, default=0.001, help="Parameter b in the gamma prior for sigma. Default is 0.001", required=False)
    parser.add_argument("--c", type=float, default=0.001, help="Parameter c for group implementation. Default is 0.001", required=False)
    parser.add_argument("--d", type=float, default=0.001, help="Parameter d for group implementation. Default is 0.001", required=False)
    parser.add_argument("--e", type=float, default=1.0, help="Parameter e for lambda. Default is 1", required=False)
    parser.add_argument("--f", type=float, default=0.5, help="Parameter f for lambda. Default is 0.5", required=False)
    parser.add_argument("--n_iter", type=int, default=1000, help="Number of maximum iterations. Default is 1000", required=False)
    parser.add_argument("--chrom", nargs='+', type=int, default=range(1,23), help="The chromosome on which the model is fitted. Default is iterating through 22 autosomes", required=False)
    parser.add_argument("--beta_std", type=str2bool, nargs='?', const=True, default=False, help="Boolean. If True, return standardized posterior SNP effect sizes. Default is False")
    parser.add_argument("--thres", type=float, default=0.001, help="Threshold value for convergence. Default is 0.001", required=False)
    parser.add_argument("--ignore_blk", type=str2bool, nargs='?', const=True, default=False, help="Boolean. If false, elements categorized as same set but in different ld will be treated as different sets. Default is False")
    parser.add_argument("--verbose", type=str2bool, nargs='?', const=True, default=False, help="Boolean. Default is False")
    parser.add_argument("--fixtau", type=str2bool, nargs='?', const=True, default=False, help="Boolean. If false, will estimate tau in algorithm, otherwise will fix at 1. Default is False")
    parser.add_argument("--flip", type=str2bool, nargs='?', const=True, default=False, help="Boolean. If false, signs of beta in summary statistics will not be flipped when there's a switch in A1 and A2 alleles. Default is False")

    param_dict = parser.parse_args()
    return param_dict


def main():
    print('='*30 + " PRSbils " +'='*30 + '\n')
    param_dict = parse_param()
    print(param_dict)
    print('-'*69 + '\n')

    for chrom in param_dict.chrom:
        print('##### process chromosome %d #####' % int(chrom))
        
        tic = time.time()
        ref_dict = parseg.parse_ref(param_dict.ref_dir + '/snpinfo_1kg_hm3', int(chrom))
        vld_dict = parseg.parse_bim(param_dict.bim_prefix, int(chrom))
        sst_dict, comm_snp = parseg.parse_sumstats(ref_dict, vld_dict, param_dict.sst_file, param_dict.n_gwas, param_dict.flip)
        set_dict, noanno_snp = parseg.parse_map(param_dict.map_file, param_dict.map_gcol, int(chrom), comm_snp)
        ld_blk_set, snp_blk_set, blk_size_set, blk_size = parseg.parse_ldblk_set(param_dict.ref_dir, sst_dict, set_dict, int(chrom), param_dict.ignore_blk, comm_snp, noanno_snp) # creating new set for snps w/o annotation
        toc = time.time()
        print('Data-parsing completed in %0.4f seconds' % (toc - tic))
        print('-'*69 + '\n')

        groups = set(set_dict['GROUP'])

        vbobj = vbclass_mcmc.VBCS(sst_dict, ld_blk_set, snp_blk_set, blk_size_set, param_dict.n_gwas, groups, param_dict.a, param_dict.b, param_dict.c, param_dict.d, param_dict.e, param_dict.f, chrom)

        tic = time.time()
        his = vbobj.fit_mcmc(param_dict.n_iter, param_dict.thres, param_dict.verbose, 500, 1.0, param_dict.fixtau)
        toc = time.time()
        print('Parameter estimation completed in %0.4f seconds' % (toc-tic))
        print('-'*69 + '\n')

        vbobj.write_beta_mcmc(param_dict.out_dir,param_dict.beta_std)
        print('Output weight file saved to\n%s' % param_dict.out_dir)
        print('='*69 + '\n')


if __name__ == '__main__':
    main()
