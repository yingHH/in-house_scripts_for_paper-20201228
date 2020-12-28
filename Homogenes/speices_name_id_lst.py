# -*- coding: utf-8 -*-
"""
Created on 2018-11-29 15:27:16
Last Modified on 2018-11-29 15:27:16

Search species id by species name and vice versa.

@Author: Ying Huang
"""
import pandas as pd


class SpNameIdLst():
    """Search species id by species name and vice versa.
    
    taxid_taxname_path: <Path> taxid_taxname contains 'taxonomy_id' and 'taxonomy_name'."""

    def __init__(self, taxid_taxname_path):
        tmp_df = pd.read_csv(taxid_taxname_path, header=None, sep='\t')
        tmp_df.columns = ['taxonomy_id', 'taxonomy_name']
        self._taxid_taxname_df = tmp_df.copy()

    def _search_by_id(self, sp_id):

        try:
            return str(self._taxid_taxname_df[
                (self._taxid_taxname_df['taxonomy_id'].astype(int) == sp_id)
            ]['taxonomy_name'].values[0])
        except:
            raise Exception('<ERR:> taxonomy_id [{}] not found.'.format(str(sp_id)))
        
    def _search_by_name(self, sp_name):
        try:
            return int(self._taxid_taxname_df[
                self._taxid_taxname_df['taxonomy_name'] == sp_name
            ]['taxonomy_id'].values[0])
        except:
            raise Exception(
                '<ERR:> taxonomy_name [{}] not found.'.format(sp_name))

    def search(self, sp_in):
        if isinstance(sp_in, str):
            return self._search_by_name(sp_in)
        elif isinstance(sp_in, int):
            return self._search_by_id(sp_in)
         
