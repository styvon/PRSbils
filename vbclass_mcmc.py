#!/usr/bin/env python3

from pyplink import PyPlink
import time
import numpy as np
import pandas as pd
import scipy as sp
import scipy.special as sc
from scipy import random
from sklearn import preprocessing
import copy
import gigrnd

class VBCS(object):
    """
    Class for bilevel continuous shrinkage
    """
    def __init__(self, sst_dict, ld_blk_set, snp_blk_set, blk_size_set, n_gwas, groups, a, b, c, d, e, f, chrom):
        self.sst_dict = sst_dict
        self.ld_blk_set = ld_blk_set
        self.snp_blk_set = snp_blk_set
        self.blk_size_set = blk_size_set
        self.n_gwas = n_gwas
        self.groups = groups
        self.n_group = len(groups)
        self.n_block = len(blk_size_set)
        # get unique names for sets
        set_names = []
        for oneblk in blk_size_set:
            set_names = []
            for oneblk in blk_size_set:
                for ss in oneblk.keys():
                    set_names.append(ss)
        self.set_names = set(set_names)
        self.n_set =  [len(oneblk.keys()) for oneblk in blk_size_set]# blk[ number of sets[]]
        self.set_size = [[[value[gg] for gg in self.groups] for key, value in oneblk.items()] for oneblk in blk_size_set]# blocks[ sets[ groups[]]]
        self.set_size_nog = [[sum([value[gg] for gg in self.groups]) for key, value in oneblk.items()] for oneblk in blk_size_set] # blocks[ sets[ ]]

        self.set_blk_map = {ss:[] for ss in self.set_names} # {set:[blk_idx]}
        self.set_size_nognob = {ss:0 for ss in self.set_names}
        self.set_size_nob = {ss:{gg:0 for gg in self.groups} for ss in self.set_names}
        for oneset in self.set_names:
            for blkid, oneblk in enumerate(blk_size_set):
                if len(self.set_size[blkid])>0:
                    if oneset in oneblk.keys():
                        self.set_blk_map[oneset].append(blkid)
                        self.set_size_nognob[oneset]+=1
                        for gg in self.groups:
                            self.set_size_nob[oneset][gg]+= oneblk[oneset][gg] # oneblk={set:{group:size}}
        self.a0 = a
        self.b0 = b
        # # for group implementation
        # self.c0 = {gg:c for gg in groups}
        # self.d0 = {gg:d for gg in groups}
        self.e0 = e # lambda params
        self.f0 = f
        self.chrom = chrom
        self.ceil = 1e30 
        self.floor =1e-30

    def fit_mcmc(self, n_iter, thres, verbose=False, burnin=500, rho=1.0, fixtau=False):
        tic = time.time()
        self._init_mcmc(verbose)

        toc = time.time()
        if verbose:
            print('Initiation completed in %0.4f seconds' % (toc-tic))
        n_iter2 = n_iter - burnin
        for iter_ in range(1, n_iter+1):
            if verbose and (iter_ % 100 == 0):
                print('==== iter %d ====' % iter_)
            # Gibbs
            tic = time.time()
            for iblk in range(self.n_block):
                self._update_beta_params_mcmc(iblk) # beta, b, delta
            temp_a = np.sum(self.a)-0.5*self.n_gwas*(len(self.a)-1)
            temp_b = max(0.5*self.n_gwas+np.sum(self.bb0)+0.5*self.n_gwas*np.sum(self.quadb), np.sum(self.bb))
            self.sig2 = 1.0/random.gamma(temp_a, 1.0/temp_b)

            for iblk in range(self.n_block):
                # update lambda2, delta, psiinv, w/ block
                self._update_lambda2_mcmc(iblk, rho = rho)
                # update posterior
                if (iter_ > burnin):
                    self.sig2_post = self.sig2_post + self.sig2/n_iter2
                    for iset,oneset in enumerate(self.set_names): 
                        for ig,gg in enumerate(self.groups):
                            if iblk in self.set_blk_map[oneset]:
                                self.beta_post[iblk][oneset][gg] = self.beta_post[iblk][oneset][gg] + self.beta[iblk][oneset][gg]/n_iter2
                                self.lambda2_post[iblk][oneset][gg] = self.lambda2_post[iblk][oneset][gg] + self.lambda2[iblk][oneset][gg]/n_iter2
            for iset,oneset in enumerate(self.set_names):
                p = self.nsnps_nognob[oneset]
                for ig,gg in enumerate(self.groups):
                    # update tau
                    if fixtau is not False:
                        self.tau2[iset][gg] = fixtau
                    else:
                        w = random.gamma(1.0, 1.0/(self.tau2[iset][gg]+1.0))
                        self.tau2[iset][gg] = random.gamma(p*self.f0+0.5, 1.0/(self.delta_sum[iset][gg]+w))
                    self.delta_sum[iset][gg] = 0.0
                    # update posterior
                    if (iter_ > burnin):
                        self.tau2_post[iset][gg] = self.tau2_post[iset][gg] + self.tau2[iset][gg]/n_iter2
            toc = time.time()
            if verbose and (iter_ % 100 == 1):
                print('Gibbs sampler: %0.4f seconds' % (toc-tic))
        print('tau2:')
        print(self.tau2_post)
        return [0, 0]

    def _init_mcmc(self, verbose=False):
        self.sig2 = self.a0/self.b0
        self.sig2_post = 0.0
        if verbose:
            print('beta...')
        self._init_betas() # blk[{set: {group:value}]}
        if verbose:
            print('psiinv...')
        self._init_psiinv()
        if verbose:
            print('quadb_YYXY...')
        self._init_quadb(verbose)
        ## group implementation
        # self.c = [{gg:0 for gg in self.groups} for kk in self.set_names]
        # self.d = [{gg:self.d0[gg] for gg in self.groups}  for oneset in self.set_names] # {set:group[]}
        self.a = [0.5*(self.n_gwas + self.nsnps_nognob[kk]) for kk in self.set_names]
        self.b = [0.0 for kk in self.set_names]

    def _init_betas(self):
        self.beta0 = []     
        self.beta = []     
        self.beta_post = [] 
        self.nsnps_nognob = {ss:0 for ss in self.set_names}    
        for iblk in range(self.n_block):
            self.beta0.append({})
            self.beta_post.append({})
            for ss in self.blk_size_set[iblk].keys():
                tempbeta = {gg:[] for gg in self.groups }
                tempbeta_zeros = {gg:[] for gg in self.groups }
                for gg in self.groups:
                    tempi = [self.sst_dict['SNP'].index(snp) for snp in self.snp_blk_set[iblk][ss][gg]]
                    tempbeta[gg] = np.atleast_1d([self.sst_dict['BETA'][ii] for ii in tempi])
                    tempbeta_zeros[gg] = np.atleast_1d([0.0 for ii in tempi])
                self.beta0[iblk].update({ss: copy.deepcopy(tempbeta)})
                self.beta_post[iblk].update({ss: copy.deepcopy(tempbeta_zeros)})
                self.nsnps_nognob[ss] += np.sum([len(betavec) for betavec in tempbeta.values()])
            tempdic = copy.deepcopy(self.beta0[iblk])
            self.beta.append(copy.deepcopy(tempdic)) # mean beta = beta*

    def _init_quadb(self, verbose=False):
        self.quadb = [0 for kk in self.set_names ]
        self.bb0 = [0.0 for kk in self.set_names] # beta*beta0
        self.bb = [0.0 for kk in self.set_names] # beta*beta

    def _init_psiinv(self):
        self.tau2 = [{gg:1.0 for gg in self.groups} for kk in self.set_names]
        self.delta_sum = [{gg:0.0 for gg in self.groups} for kk in self.set_names] # for tau2 update
        self.tau2_post = [{gg:0.0 for gg in self.groups} for kk in self.set_names]
        self.lambda2 = []
        self.lambda2_post = []
        self.psiinv = []
        for iblk in range(self.n_block):
            self.lambda2.append({})
            self.lambda2_post.append({})
            self.psiinv.append({})
            for ss in self.blk_size_set[iblk].keys():
                templambda2 = {gg:[] for gg in self.groups }
                templambda2_zeros = {gg:[] for gg in self.groups }
                temppsiinv = {gg:[] for gg in self.groups }
                for gg in self.groups:
                    tempi = [self.sst_dict['SNP'].index(snp) for snp in self.snp_blk_set[iblk][ss][gg]]
                    templambda2[gg] = np.atleast_1d([1.0 for ii in tempi])
                    templambda2_zeros[gg] = np.atleast_1d([0.0 for ii in tempi])
                    temppsiinv[gg] = np.atleast_1d([1.0 for ii in tempi])
                self.lambda2[iblk].update({ss:templambda2})
                self.lambda2_post[iblk].update({ss:templambda2_zeros})
                self.psiinv[iblk].update({ss:temppsiinv})

    def _update_beta_params_mcmc(self,iblk): 
        for iset, oneset in enumerate(self.set_names):
            sig2 = self.sig2
            if oneset in self.blk_size_set[iblk].keys():
                if iblk == self.set_blk_map[oneset][0]:
                    self.quadb[iset] = 0.0
                    self.bb0[iset] = 0.0
                    self.bb[iset] = 0.0
                    self.b[iset] = 0.5*self.n_gwas
                for gg in self.groups:
                    LD = self.ld_blk_set[iblk][oneset][gg]
                    var_pre = LD+sp.diag(1.0/self.lambda2[iblk][oneset][gg])
                    var_pre2 = sp.linalg.cholesky(var_pre) #dinvt_chol
                    tempbeta = sp.linalg.solve_triangular(var_pre2, self.beta0[iblk][oneset][gg], trans='T') + sp.sqrt(sig2/self.n_gwas)*random.randn(len(self.lambda2[iblk][oneset][gg]))
                    self.beta[iblk][oneset][gg] = np.atleast_1d(sp.linalg.solve_triangular(var_pre2, tempbeta, trans='N'))
                    tempquadb = sp.dot(sp.dot(self.beta[iblk][oneset][gg].T, var_pre), self.beta[iblk][oneset][gg])
                    self.quadb[iset] += tempquadb
                    self.bb0[iset] += self.n_gwas*(0.0-np.sum(self.beta[iblk][oneset][gg]*self.beta0[iblk][oneset][gg]))
                    self.bb[iset] += self.n_gwas*0.5*np.sum(self.beta[iblk][oneset][gg]**2/self.lambda2[iblk][oneset][gg])
    
    def _update_lambda2_mcmc(self,iblk,rho=1.0): 
        for iset, oneset in enumerate(self.set_names):
            sig2 = self.sig2
            if oneset in self.blk_size_set[iblk].keys():
                p = self.nsnps_nognob[oneset]
                for ig,gg in enumerate(self.groups):
                    p_blk = len(self.beta0[iblk][oneset][gg])
                    # delta 
                    delta = random.gamma(self.e0+self.f0, 1.0/(self.lambda2[iblk][oneset][gg]+self.tau2[iset][gg]))
                    self.delta_sum[iset][gg] += np.sum(delta)
                    # lambda2
                    for j in range(p_blk):
                        self.lambda2[iblk][oneset][gg][j] = gigrnd.gigrnd(self.e0-0.5, 2.0*delta[j], self.n_gwas*self.beta[iblk][oneset][gg][j]**2/self.sig2)
                    self.lambda2[iblk][oneset][gg][self.lambda2[iblk][oneset][gg]>rho] = rho
    
    def write_beta_mcmc(self, out_dir, beta_std):
        out_file = out_dir + '_beta_chr%s.txt' % (self.chrom)
        with open(out_file, 'w') as ff:
            for iblk in range(self.n_block):
                for oneset in self.beta_post[iblk].keys():
                    for gg in self.beta_post[iblk][oneset].keys():
                        sst_id = [self.sst_dict['SNP'].index(snp) for snp in self.snp_blk_set[iblk][oneset][gg]]
                        for ii,sid in enumerate(sst_id):
                            chrom = self.chrom
                            snp = self.sst_dict['SNP'][sid]
                            bp = self.sst_dict['BP'][sid]
                            a1 = self.sst_dict['A1'][sid]
                            a2 = self.sst_dict['A2'][sid]
                            maf = self.sst_dict['MAF'][sid]
                            beta0 = self.beta0[iblk][oneset][gg][ii]
                            beta = self.beta_post[iblk][oneset][gg][ii]
                            # convert beta
                            if beta_std == 'False':
                                beta /= sp.sqrt(2.0*maf*(1.0-maf))
                            ff.write('%s\t%s\t%d\t%s\t%s\t%.6e\t%.6e\t%s\t%s\t%.4f\n' % (chrom, snp, bp, a1, a2, beta0, beta, oneset, gg, maf))
    


