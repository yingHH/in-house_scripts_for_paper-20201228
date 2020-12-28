# -*- coding: utf-8 -*-
"""
Created on 2018-11-29 09:39:13
Last Modified on 2018-11-29 09:39:13

Find homologous genes based on './homo_db' or other customer database.
(database download from: ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current).

@Author: Ying Huang
"""
import os
import pandas as pd


def find_homogenes(homo_db, species_ids):
    """Find homologous genes based on './homo_db' or other customer database.
    (database download from: https://ftp.ncbi.nih.gov/pub/HomoloGene/current/).
    
    homo_db: <Path> database directory contained files 'homologene.data, 
        taxid_taxname, all_proteins.data'. All of the file can be download
        from 'https://ftp.ncbi.nih.gov/pub/HomoloGene/current/'.
    species_ids: <iter> a list of species ids."""

    homologene_data_path = os.path.join(homo_db, 'homologene.data')
    homo_data = pd.read_csv(homologene_data_path, header=None, sep='\t')
    homo_data.columns = ['group_id', 'taxonomy_id', 'gene_id', 'gene_name', 'protein_gi', 'protein_id']

    # get each species data by species id
    species_dfs = []
    for sp_id in species_ids:
        sp_id = int(sp_id)
        # get genes list of a species id ('sp_id')
        tmp_df = homo_data[homo_data['taxonomy_id'].astype(int) == sp_id]
        # add species id to each columns 
        tmp_df.columns = ['{}_{}'.format(n, str(sp_id)) for n in tmp_df.columns]

        # only save 'group_id', 'gene_id' and 'protein_id' columns
        species_dfs.append(
            tmp_df.iloc[:, [0, 3, 5]].copy())

    # merge all species data
    mreged_data = pd.DataFrame()
    for i, df in enumerate(species_dfs):
        if i == 0:
            mreged_data = df
            continue
        # merge df based on group id
        mreged_data = pd.merge(
            left=mreged_data, left_on=mreged_data.columns[0], 
            right=df, right_on=df.columns[0], 
            how='outer',
        )

    # remove umatched rows
    na_mask = mreged_data.apply(lambda x: x.isna().any(), axis=1)
    mreged_data = mreged_data[~(na_mask)]

    return mreged_data.copy()

    




