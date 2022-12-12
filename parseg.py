import sys
import scipy as sp
from scipy.stats import norm
import h5py
import gc
from operator import itemgetter
import utils # mergeDict
from metax import PredictionModel

P_MIN = 1E-323

def parse_map(map_file, map_gcol, chrom, include_snp = None):
    """
    Parse annotation map file.

    :param map_file: annotation map file
    :param map_gcol: column id for group information, integer or None
    :param chrom: chromosome to be included
    :param include_snp: list of rsids to be included. Default is None (include all snps)
    """
    print('... parse map file: %s ...' % (map_file))
    set_dict = {'SNP':[], 'SET':[], 'GROUP':[]}
    anno_snp =set([])
    if include_snp is not None:
        # include only specified snps
        if map_gcol== None:    
            # if group column unspecified, use default '0'
            with open(map_file) as ff:
                for line in ff:
                    ll = (line.strip()).split()
                    if (int(ll[0]) == chrom) & (ll[1] in include_snp):
                        set_dict['SNP'].append(ll[1])
                        set_dict['SET'].append(ll[2])
                        set_dict['GROUP'].append('0')
                        anno_snp.add(ll[1])
        else:
            with open(map_file) as ff:
                for line in ff:
                    ll = (line.strip()).split()
                    if (int(ll[0]) == chrom) & (ll[1] in include_snp):
                        set_dict['SNP'].append(ll[1])
                        set_dict['SET'].append(ll[2])
                        set_dict['GROUP'].append(ll[map_gcol-1])
                        anno_snp.add(ll[1])
    else:
        # include all snps in map file
        if map_gcol== None:    
            # if group column unspecified, use default '0'
            with open(map_file) as ff:
                for line in ff:
                    ll = (line.strip()).split()
                    if int(ll[0]) == chrom:
                        set_dict['SNP'].append(ll[1])
                        set_dict['SET'].append(ll[2])
                        set_dict['GROUP'].append('0')
        else:
            with open(map_file) as ff:
                for line in ff:
                    ll = (line.strip()).split()
                    if int(ll[0]) == chrom:
                        set_dict['SNP'].append(ll[1])
                        set_dict['SET'].append(ll[2])
                        set_dict['GROUP'].append(ll[map_gcol-1])
        include_snp = set(set_dict['SNP'])
    noanno_snp = include_snp - anno_snp
    if len(noanno_snp)>0:
        for temp_snp in noanno_snp:
            set_dict['SNP'].append(temp_snp)
            set_dict['SET'].append('NO_ANNO')
            set_dict['GROUP'].append('0')
    else:
        noanno_snp = set([])

    print('... %d SNPs on chromosome %d read from %s ...' % (len(set_dict['SNP']), chrom, map_file))
    return set_dict, noanno_snp


def parse_ldblk_set(ldblk_dir, sst_dict, set_dict, chrom, ignore_blk=False, comm_snp=None, noanno_snp=None):
    """
    Parse LD block info by annotation sets.

    :param ldblk_dir: LD block data directory
    :param sst_dict: dictionary parsed from summary statistics
    :param set_dict: dictionary parsed from annotation map
    :param chrom: chromosome to be included
    :param ignore_blk: whether grouped by block. Default is False
    :param comm_snp: whether grouped by block. Default is False
    :param noanno_snp: if noanno_snp=None, variants w/o annotation will be dropped
    """

    if comm_snp is None:
        comm_snp = [a for a in sst_dict['SNP'] if a in set_dict['SNP']]
    print('... %d common SNPs in the reference, sumstats, validation and set info set ...' % len(comm_snp))
    if noanno_snp is None:
        noanno_snp = set(sst_dict['SNP'])-set(set_dict['SNP'])
    print('... %d SNPs not annotated ...' % len(noanno_snp))

    print('... parse reference LD with set info on chromosome %d ...' % chrom)

    chr_name = ldblk_dir + '/ldblk_1kg_chr' + str(chrom) + '.hdf5'
    hdf_chr = h5py.File(chr_name, 'r')
    n_blk = len(hdf_chr)
    ld_blk = [sp.array(hdf_chr['blk_'+str(blk)]['ldblk']) for blk in range(1,n_blk+1)]
    if sys.version_info.major == 2:
        snp_blk = [list(hdf_chr['blk_'+str(blk)]['snplist']) for blk in range(1,n_blk+1)]
    elif sys.version_info.major ==3:
        snp_blk = []
        for blk in range(1,n_blk+1):
            temp = list(hdf_chr['blk_'+str(blk)]['snplist'])
            snp_blk.append([onetemp.decode("utf-8") for onetemp in temp])
    
    # initialize dict{set: snp}
    sets = set(set_dict['SET'])
    n_sets = len(sets)
    groups = set(set_dict['GROUP'])
    n_grp = len(groups)
    # set: group: snp, unequal length per item
    set_dict2 = {new_list: {gg:[] for gg in groups} for new_list in sets} 

    for (ii, snp) in enumerate(set_dict['SNP']):
        gg = set_dict['GROUP'][ii]
        if snp in comm_snp:
            oneset = set_dict['SET'][ii]
            set_dict2[oneset][gg].append(snp)

    # initialize outputs: block[ set{ group{}}]
    ld_blk_set = [{new_set: {gg: [] for gg in groups}  for new_set in sets} for blk in range(n_blk)] # list of n_blk dicts
    snp_blk_set = [{new_set: {gg: [] for gg in groups} for new_set in sets} for blk in range(n_blk)] # list of n_blk dicts
    blk_size_set = [{new_set: {gg: 0 for gg in groups} for new_set in sets}  for blk in range(n_blk)] # list of n_blk dicts
    blk_size = [] # inx for non-empty blocks

    # for each set within each ld block, get ld, snps, size; get accum size per block  
    for blk in range(n_blk):
        oneblk_size = 0
        for oneset in sets:
            oneset_size = []
            idx = {gg: [] for gg in groups} # snp index in ld ref
            idx_sst = {gg: [] for gg in groups} # snp index in sst
            for (gi, onegroup) in enumerate(groups):
                onegroup_size = 0
                for (ii, snp) in enumerate(snp_blk[blk]):                      
                    if snp in set_dict2[oneset][onegroup]:
                        idx[onegroup].append(ii)
                        onegroup_size += 1
                        idx_sst[onegroup].append(sst_dict['SNP'].index(snp))
                if onegroup_size == 0:
                    temp = ld_blk_set[blk][oneset].pop(onegroup, None) # remove empty set from the blk dictionary    
                    temp = snp_blk_set[blk][oneset].pop(onegroup, None) # remove empty set from the blk dictionary    
                    temp = blk_size_set[blk][oneset].pop(onegroup, None) # remove empty set from the blk dictionary    
                elif onegroup_size == 1:
                    # ld per set
                    ld_blk_set[blk][oneset][onegroup]=ld_blk[blk][sp.ix_(idx[onegroup],idx[onegroup])]
                    # snps in one set
                    snp_blk_set[blk][oneset][onegroup]=[snp_blk[blk][idx[onegroup][0]]]
                    # size of set
                    blk_size_set[blk][oneset][onegroup] = onegroup_size
                else:
                    # ld per set
                    flip = [sst_dict['FLP'][jj] for jj in idx_sst[onegroup]]
                    ld_blk_set[blk][oneset][onegroup]=ld_blk[blk][sp.ix_(idx[onegroup],idx[onegroup])]*sp.outer(flip,flip)
                    _, s, v = sp.linalg.svd(ld_blk_set[blk][oneset][onegroup])
                    h = sp.dot(v.T, sp.dot(sp.diag(s), v))
                    ld_blk_set[blk][oneset][onegroup] = (ld_blk_set[blk][oneset][onegroup]+h)/2           
                    # snps in one set
                    snp_blk_set[blk][oneset][onegroup]=list(itemgetter(*idx[onegroup])(snp_blk[blk]))
                    # size of set
                    blk_size_set[blk][oneset][onegroup] = onegroup_size
                oneset_size.append(onegroup_size)
            if sum(oneset_size) == 0:
                temp = ld_blk_set[blk].pop(oneset, None) # remove empty set from the blk dictionary    
                temp = snp_blk_set[blk].pop(oneset, None) # remove empty set from the blk dictionary    
                temp = blk_size_set[blk].pop(oneset, None) # remove empty set from the blk dictionary    
            oneblk_size += sum(oneset_size)
        blk_size.append(oneblk_size)

    # if ignore ld block info, combine sets with same name across blocks
    if ignore_blk:
        ld_blk_set_ignoreblk = ld_blk_set[0]
        snp_blk_set_ignoreblk = snp_blk_set[0]
        blk_size_set_ignoreblk = blk_size_set[0]
        blk_size_ignoreblk = None

        for blk in range(1,n_blk):
            if blk_size[blk]!=0:
                ld_blk_set_ignoreblk = utils.mergeDict(ld_blk_set_ignoreblk, ld_blk_set[blk], mode="matrix", groups = groups)
                snp_blk_set_ignoreblk = utils.mergeDict(snp_blk_set_ignoreblk, snp_blk_set[blk], mode="list", groups = groups)
                blk_size_set_ignoreblk = utils.mergeDict(blk_size_set_ignoreblk, blk_size_set[blk], mode="numeric", groups = groups)
        return ld_blk_set_ignoreblk, snp_blk_set_ignoreblk, blk_size_set_ignoreblk, blk_size_ignoreblk

    return ld_blk_set, snp_blk_set, blk_size_set, blk_size


def parse_ref(ref_file, chrom):
    """ Title: PRScs source code, parse_ref
    Author: Tian Ge
    Date: 2021
    Code version: 1.0.0
    Availability: https://github.com/getian107/PRScs """

    print('... parse reference file: %s ...' % ref_file)

    ref_dict = {'CHR':[], 'SNP':[], 'BP':[], 'A1':[], 'A2':[], 'MAF':[]}
    with open(ref_file) as ff:
        header = next(ff)
        for line in ff:
            ll = (line.strip()).split()
            if int(ll[0]) == chrom:
                ref_dict['CHR'].append(chrom)
                ref_dict['SNP'].append(ll[1])
                ref_dict['BP'].append(int(ll[2]))
                ref_dict['A1'].append(ll[3])
                ref_dict['A2'].append(ll[4])
                ref_dict['MAF'].append(float(ll[5]))

    print('... %d SNPs on chromosome %d read from %s ...' % (len(ref_dict['SNP']), chrom, ref_file))
    return ref_dict


def parse_bim(bim_file, chrom):
    """ Title: PRScs source code, parse_bim
    Author: Tian Ge
    Date: 2021
    Code version: 1.0.0
    Availability: https://github.com/getian107/PRScs """

    print('... parse bim file: %s ...' % (bim_file + '.bim'))

    vld_dict = {'SNP':[], 'A1':[], 'A2':[]}
    with open(bim_file + '.bim') as ff:
        for line in ff:
            ll = (line.strip()).split()
            if int(ll[0]) == chrom:
                vld_dict['SNP'].append(ll[1])
                vld_dict['A1'].append(ll[4])
                vld_dict['A2'].append(ll[5])

    print('... %d SNPs on chromosome %d read from %s ...' % (len(vld_dict['SNP']), chrom, bim_file + '.bim'))
    return vld_dict

def parse_sumstats(ref_dict, vld_dict, sst_file, n_subj, flip=True):
    """
    Parse summary statistics file.
    (Edited based on PRScs source code)
    """

    print('... parse sumstats file: %s ...' % sst_file)
    sst_dict = {'SNP':[], 'A1':[], 'A2':[]}
    # check sst_file (summary statistics file) format
    sst_format = "default"
    if sst_file.endswith(".db"):
        sst_format = "db"

    if sst_format == "default":
        # default plain text format
        with open(sst_file) as ff:
            header = next(ff)
            for line in ff:
                ll = (line.strip()).split()
                sst_dict['SNP'].append(ll[0])
                sst_dict['A1'].append(ll[1])
                sst_dict['A2'].append(ll[2])
    elif sst_format == "db":
        # predictdb format
        sst_model = PredictionModel.load_model(sst_file)

    print('... %d SNPs read from %s ...' % (len(sst_dict['SNP']), sst_file))

    mapping = {'A': 't', 'T': 'a', 'C': 'g', 'G': 'c'}
    ref_snp = set(zip(ref_dict['SNP'], ref_dict['A1'], ref_dict['A2'])) | \
              set(zip(ref_dict['SNP'], ref_dict['A2'], ref_dict['A1'])) | \
              set(zip(ref_dict['SNP'], utils.multiReplace(ref_dict['A1'], mapping), utils.multiReplace(ref_dict['A2'], mapping) )) | \
              set(zip(ref_dict['SNP'], utils.multiReplace(ref_dict['A2'], mapping), utils.multiReplace(ref_dict['A1'], mapping) ))

    vld_snp = set(zip(vld_dict['SNP'], vld_dict['A1'], vld_dict['A2']))

    sst_snp = set(zip(sst_dict['SNP'], sst_dict['A1'], sst_dict['A2'])) | \
              set(zip(sst_dict['SNP'], sst_dict['A2'], sst_dict['A1'])) | \
              set(zip(sst_dict['SNP'], utils.multiReplace(ref_dict['A1'], mapping), utils.multiReplace(ref_dict['A2'], mapping))) | \
              set(zip(sst_dict['SNP'], utils.multiReplace(ref_dict['A2'], mapping), utils.multiReplace(ref_dict['A1'], mapping)))

    comm_snp = ref_snp & vld_snp & sst_snp

    print('... %d common SNPs in the reference, sumstats, and validation set ...' % len(comm_snp))

    n_sqrt = sp.sqrt(n_subj)
    sst_eff = {}
    with open(sst_file) as ff:
        header = (next(ff).strip()).split()
        header = [col.upper() for col in header]
        for line in ff:
            ll = (line.strip()).split()
            snp = ll[0]; a1 = ll[1]; a2 = ll[2]
            if flip:
                if ((snp, a1, a2) in comm_snp) or ((snp, utils.multiReplace1(a1, mapping), utils.multiReplace1(a2, mapping)) in comm_snp):
                    if 'BETA' in header:
                        beta = float(ll[3])
                    elif 'OR' in header:
                        beta = sp.log(float(ll[3]))

                    p = max(float(ll[4]), P_MIN)
                    beta_std = sp.sign(beta)*abs(norm.ppf(p/2.0))/n_sqrt
                    sst_eff.update({snp: beta_std})
                elif ((snp, a2, a1) in comm_snp) or ((snp, utils.multiReplace1(a2, mapping), utils.multiReplace1(a1, mapping)) in comm_snp):
                    if 'BETA' in header:
                        beta = float(ll[3])
                    elif 'OR' in header:
                        beta = sp.log(float(ll[3]))

                    p = max(float(ll[4]), P_MIN)
                    beta_std = -1*sp.sign(beta)*abs(norm.ppf(p/2.0))/n_sqrt
                    sst_eff.update({snp: beta_std})
            else:
                if 'BETA' in header:
                    beta = float(ll[3])
                elif 'OR' in header:
                    beta = sp.log(float(ll[3]))

                p = max(float(ll[4]), P_MIN)
                beta_std = sp.sign(beta)*abs(norm.ppf(p/2.0))/n_sqrt
                sst_eff.update({snp: beta_std})
            # end if flip

    sst_dict = {'CHR':[], 'SNP':[], 'BP':[], 'A1':[], 'A2':[], 'MAF':[], 'BETA':[], 'FLP':[]}
    for (ii, snp) in enumerate(ref_dict['SNP']):
        if snp in sst_eff:
            sst_dict['SNP'].append(snp)
            sst_dict['CHR'].append(ref_dict['CHR'][ii])
            sst_dict['BP'].append(ref_dict['BP'][ii])
            a1 = ref_dict['A1'][ii]
            a2 = ref_dict['A2'][ii]

            if flip:
                if (snp, a1, a2) in comm_snp:
                    sst_dict['A1'].append(a1)
                    sst_dict['A2'].append(a2)
                    sst_dict['MAF'].append(ref_dict['MAF'][ii])
                    sst_dict['FLP'].append(1)
                elif (snp, a2, a1) in comm_snp:
                    sst_dict['A1'].append(a2)
                    sst_dict['A2'].append(a1)
                    sst_dict['MAF'].append(1-ref_dict['MAF'][ii])
                    sst_dict['FLP'].append(-1)
                elif (snp, utils.multiReplace1(a1, mapping), utils.multiReplace1(a2, mapping)) in comm_snp:
                    sst_dict['A1'].append(utils.multiReplace1(a1, mapping))
                    sst_dict['A2'].append(utils.multiReplace1(a2, mapping))
                    sst_dict['MAF'].append(ref_dict['MAF'][ii])
                    sst_dict['FLP'].append(1)
                elif (snp, utils.multiReplace1(a2, mapping), utils.multiReplace1(a1, mapping)) in comm_snp:
                    sst_dict['A1'].append(utils.multiReplace1(a2, mapping))
                    sst_dict['A2'].append(utils.multiReplace1(a1, mapping))
                    sst_dict['MAF'].append(1-ref_dict['MAF'][ii])
                    sst_dict['FLP'].append(-1)
            else:
                sst_dict['A1'].append(a1)
                sst_dict['A2'].append(a2)
                sst_dict['MAF'].append(ref_dict['MAF'][ii])
                sst_dict['FLP'].append(1)
            # end if flip
            sst_dict['BETA'].append(sst_eff[snp])

    comm_snp_idonly = set([a1 for (a1,a2,a3) in comm_snp])

    return sst_dict, comm_snp_idonly