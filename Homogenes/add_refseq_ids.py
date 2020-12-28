# -*- coding: utf-8 -*-
"""
Created on 2018-11-29 14:14:33
Last Modified on 2018-11-29 14:14:33

Add refseq id list to homogenes list base on all_proteins.data file.

@Author: Ying Huang
"""
import os
import pandas as pd


def get_protein_refseq_lst_of_species(taxonomy_id, all_proteins_data):

    species_lst = all_proteins_data[
        all_proteins_data['taxonomy_id'] == taxonomy_id
    ]

    species_lst.columns = [
        '{}_{}'.format(n, taxonomy_id) for n in species_lst.columns
    ]

    return species_lst.iloc[:, [1,2]].copy()


def add_gene_ids(homo_db, homogenes_df, species_ids):
    """Add refseq id list to homogenes list base on all_proteins.data file.
    
    homo_db: <Path> database directory contained files 'homologene.data, 
        taxid_taxname, all_proteins.data'. All of the file can be download
        from 'ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current'.
    homogenes_df: <DataFrame> result of 'find_homogenes.find_homogenes(homo_db,
        species_ids)', which contains 'group_id', 'gene_name', 'protein_id' of
        homogenes.
    species_ids: <iter> a list of species ids."""

    homogenes_df = homogenes_df.copy()

    all_proteins_data_path = os.path.join(homo_db, 'all_proteins.data')
    all_proteins_data = pd.read_csv(all_proteins_data_path, header=None, sep='\t')[[0, 4, 5]]
    all_proteins_data.columns = ['taxonomy_id', 'protein_id', 'refseq_id']

    # merge refseq id to homogenes_df
    for sp_id in species_ids:
        tmp_protein_refseq_lst = get_protein_refseq_lst_of_species(
            sp_id, all_proteins_data
        )

        homogenes_df = pd.merge(
            homogenes_df, tmp_protein_refseq_lst,
            on=tmp_protein_refseq_lst.columns[0],
            how='outer',
        )

    # remove umatched rows
    na_mask = homogenes_df.apply(lambda x: x.isna().any(), axis=1)
    homogenes_df = homogenes_df[~(na_mask)]

    # sort columns by 'taxonomy_id'
    sorted_col = homogenes_df.columns.tolist()
    sorted_col.sort(
        key=lambda x: int(x.split('_')[-1])
    )

    return homogenes_df[sorted_col].copy()
