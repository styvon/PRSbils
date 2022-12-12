#!/bin/sh

pip3 install -r requirements.txt

python3 PRSbils.py --ref_dir=./1000G/ldblk_1kg_eur --bim_prefix=./data/genotypes_plink_chr1_train --sst_file=./data/sumstat.txt --map_file=./data/snpmap.txt --n_gwas=50000 --out_dir=./data/output_trained_weights_prsbils --chrom=1 --n_iter=1000 --thres=1e-04 --beta_std=False --flip=True
