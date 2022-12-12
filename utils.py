#!/usr/bin/env python3

"""
Utility functions for general data transformation.

"""
import numpy as np
import scipy.special as sc
from scipy.linalg import block_diag
import pandas as pd
import statistics
import math
from sklearn.linear_model import LogisticRegression
from itertools import combinations

def mergeDict(dict1, dict2, mode="matrix", groups=None):
    ''' 
    Merge dictionaries in order dict1, dict2
    @params    dict1
    @params dict2
    @params    mode     "matrix": form block diagonal matrix; "list": append list in dict2 to that in dict1 with same key; "numeric": sum numeric elements with same key
    @example
    test1 = {'a':[np.array([[1,2],[3,4]])]}
    test2 = {'a':[np.array([[5,6],[7,8]])]}
    test3 = utils.mergeDict(test1, test2, mode="matrix")

    test1 = {'a':[1,2], 'b':[3,4]}
    test2 = {'a':[5,6],'c':[7,8]}
    test3 = utils.mergeDict(test1, test2, mode="list")

    test1 = {'a':{'a':[1,2], 'b':[3,4]}, 'b':{'a':[1,2], 'b':[3,4]}}
    test2 = {'a':{'a':[5,6], 'b':[7,8]},'c':{'a':[9,10], 'b':[11,12]}}
    test3 = utils.mergeDict(test1, test2, mode="list", groups = set(['a','b']))
    '''
    # dict3 = {**dict1, **dict2}
    dict3 = dict2.copy()
    dict3.update(dict1)

    if groups == None:
        if mode=="matrix":
            for key, value in dict3.items():
                if key in dict1 and key in dict2:
                    dict3[key]=block_diag(dict3[key],dict2[key])
        elif mode =="list":
            for key, value in dict3.items():
                if key in dict1 and key in dict2:
                    dict3[key]= np.append(dict3[key], dict2[key])
                    # dict3[key].append(dict2[key])
        elif mode =="numeric":
            for key, value in dict3.items():
                if key in dict1 and key in dict2:
                    dict3[key] = dict3[key] + dict2[key]
    else: 
        # if groups not none
        for gg in groups:
            if mode=="matrix":
                for key, value in dict3.items():
                    if key in dict1 and key in dict2:
                        dict3[key][gg]=block_diag(dict3[key][gg],dict2[key][gg])
            elif mode =="list":
                for key, value in dict3.items():
                    if key in dict1 and key in dict2:
                        dict3[key][gg]= np.append(dict3[key][gg], dict2[key][gg])
                        # dict3[key].append(dict2[key])
            elif mode =="numeric":
                for key, value in dict3.items():
                    if key in dict1 and key in dict2:
                        dict3[key][gg] = dict3[key][gg] + dict2[key][gg]

    return dict3

def getOR(x,y, by_decile=10):
    x_array = np.squeeze(np.asarray(x))
    y_array = np.squeeze(np.asarray(y))
    if by_decile:
        deciles = pd.qcut(x_array,by_decile)
        deciles_dum = pd.get_dummies(data=deciles, drop_first=True)
        # dat = np.column_stack((x,deciles_dum))
        onevec = np.ones(len(y_array))
        dumvec = deciles_dum.values[:,-1]
        # dat = deciles_dum
        dat = np.column_stack((onevec, dumvec))
        glm_model = LogisticRegression(random_state=0, solver='lbfgs',fit_intercept=False).fit(dat, y_array)
        coef = glm_model.coef_[0]
        # out = [math.exp(c) for c in coef]
        out = math.exp(coef[-1])
        return(out)

    else:
        print('Need to specify number of deciles')
        return(0)

def getR2(x,y, by_decile=0):
    x_array = np.squeeze(np.asarray(x))
    # x_array = (x_array - mean(x_array)) / std(x_array)
    y_array = np.squeeze(np.asarray(y))
    if by_decile:
        deciles = pd.qcut(x_array,by_decile)
        cat = deciles.categories
        r2 = []
        for c in cat:
            idx = np.where(deciles==c)
            corr_matrix = np.corrcoef(x_array[idx], y_array[idx])
            corr_xy = corr_matrix[0,1]
            r2.append(corr_xy**2)
        return(r2)
    else:
        correlation_matrix = np.corrcoef(x_array, y_array)
        correlation_xy = correlation_matrix[0,1]
        r2 = correlation_xy**2
        return(r2)

def getIOU(x,y):
    """
    Intersection over Union for two sets
    """
    set1 = set(x)
    set2 = set(y)
    m1 = len(set1)
    m2 = len(set2)
    m12 = len(set1 & set2)
    iou = m12/(m1+m2-m12)
    return(iou)

def pairwise(l):
    """Get an iterable pair from list"""
    out = list(combinations(l, 2))
    return out

def _lnexp1_taylor(x, n_term):
    """Series expansion for E1"""
    r = 0; s = -1
    for i in range(n_term, 0, -1):
        r = (r + s) * i / x
        s = -s
    return r*s

def lnexp1(x, n_term=18):
    """Calculation of log of Exponential integral E1 using series expansion"""
    x = np.array(x)
    use_approx = x >= 50

    ax = np.where(use_approx, x, 100)
    sx = np.where(use_approx, 1, x)
    approx = 0.0-x -np.log(x) + np.log1p(_lnexp1_taylor(ax, n_term))
    temp = np.log(sc.exp1(sx))
    return np.where(use_approx, approx, temp) * 1

def safeexp(x, floor=-710,ceil=710):
    if x <= floor:
        out = 0
    elif x >=ceil:
        out = np.exp(ceil-1)
    else:
        out = np.exp(x)
    return out

def multiReplace(indata, dictionary):
    for key in dictionary.keys():
        indata = [aa.replace(key, dictionary[key]) for aa in indata]
    output = [aa.upper() for aa in indata]
    return output

def multiReplace1(indata, dictionary):
    for key in dictionary.keys():
        indata = indata.replace(key, dictionary[key])
    output = indata.upper()
    return output