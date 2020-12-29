# -*- coding: utf-8 -*-
"""
Author: Ying Huang

Date: 2020-05-21 16:50:17
Last Modified by: Ying Huang
Last Modified time: 2020-05-21 16:50:17

Description: 
"""
import pandas as pd
import numpy as np

## sub func

# JSD
def I(p):
    """
    p: is a gene expression probibilaty.
    """
    if p == 0:
        return 0
    else:
        return - p * np.log2(p)
    

def H(se):
    """
    se: is a gene expression distribution. It is a pd.Serise across all tisues. 
    """
    se = se[(se > 0)] # remove zeroes
    log_se = np.log2(se)
    return (- se.values * log_se.values).sum()


def jsd_spec(df):
    '''
    df: A DataFrame of all genes expression value (normalized; rows) across all tissues (columns).
    '''
    H_diff_values = df.apply(lambda x: H(x/2) - 1/2 * H(x), axis=1)
    #print('H(e/2) - 1/2 H(e):', H_diff_values)
    I_diff_values = df.applymap(lambda x: I((x + 1)/2) - I(x/2))
    #print('I diff:', I_diff_values)
    return (I_diff_values.T + H_diff_values).T
    

# Main
#       
def jsd(df):
    '''
    df: A DataFrame of all genes expression value (normalized; rows) across all tissues (columns).
    '''
    assert isinstance(df, pd.DataFrame)

    tsv_df = 1 - jsd_spec(df) ** (1/2)
    return tsv_df