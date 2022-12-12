
# PRSbils  

<!-- badges: start -->

<!-- badges: end -->

## Updates  

![Update log](NEWS.md)  

## Introduction  

Polygenic risk score with bilevel continuous shrinkage for incorporating functional annotations.  

## Usage

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
  --map_gcol MAP_GCOL   Column id for group info in map_file. Default is None
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
