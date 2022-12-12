
# PRSbils  



## Updates  

![Update log](NEWS.md)  

## Introduction  

Polygenic risk score with bilevel continuous shrinkage for incorporating functional annotations.  



## How to use

1. Download/Clone PRSbils (this repository) to your local machine
2. Install required python packages: `pip3 install -r requirements.txt`
3. Get the input files ready
  - Download and extract LD reference panel from 1000 Genomes Project phase 3 (the same procedure as in [PRSCS](https://github.com/getian107/PRScs/blob/master/README.md)):
    - [African](https://www.dropbox.com/s/mq94h1q9uuhun1h/ldblk_1kg_afr.tar.gz?dl=0 "African")
    - [American](https://www.dropbox.com/s/uv5ydr4uv528lca/ldblk_1kg_amr.tar.gz?dl=0 "American")
    - [East Asian](https://www.dropbox.com/s/7ek4lwwf2b7f749/ldblk_1kg_eas.tar.gz?dl=0 "East Asian")
    - [European](https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz?dl=0 "European")
    - [South-eastern Asian]https://www.dropbox.com/s/hsm0qwgyixswdcv/ldblk_1kg_sas.tar.gz?dl=0 "South-eastern Asian"
  - Genotype file
  - GWAS summary statistics
  - Mapping file from SNPs to annotation groups
4. Run `PRSbils.py` to get the trained weights for SNPs with different functional annotations.

## Sample command

You can test the program using the following bash command after downloading the 1000G European LD referenfce panel file under the home directory:

```
pip3 install -r requirements.txt

python3 PRSbils.py \
        --ref_dir=./1000G/ldblk_1kg_eur \
        --bim_prefix=./data/genotypes_plink_chr1_train \
        --sst_file=./data/sumstat.txt \
        --map_file=./data/snpmap.txt \
        --n_gwas=50000 \
        --out_dir=./data/output_trained_weights_prsbils \
        --chrom=1 \
        --n_iter=1000 \
        --thres=1e-04 \
        --beta_std=False \
        --flip=True

```

The output file is  `./data/output_trained_weights_prsbils_beta_chr1.txt`, which is a tab-delimited text file with the following columns:
```
[Chromosome] [SNP ID] [SNP Position] [Allele 1] [Allele 2] [Weight in the summary statistics file] [Weight from PRSbils] [Annotation Group] [Subgroup (not yet implemented)] [Allele Frequency]
```

Please check the next section or run `python3 PRSbils.py --help` for details about the parameters.

## Details about the parameters

```
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
```
