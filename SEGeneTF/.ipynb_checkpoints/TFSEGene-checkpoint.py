# -*- coding: utf-8 -*-
"""
Created on 2018-12-01 01:48:37
Last Modified on 2018-12-01 01:48:37

Read in TF file of fimo, and assign SEs to genes based on it.

@Author: Ying Huang
"""
from itertools import chain
import os
import pandas as pd


WorkDir = os.path.curdir
PWD = os.path.split(os.path.realpath(__file__))[0]
TF_motif_dict_path = os.path.join(PWD, 'motif_db', 'MotifDictionary.txt')


def count_bind_freq_se(fimo_path, se_gene, freq_thread=3):
    """Read in fimo output (*.fimo.txt), then count binding frequence of the TF.
    
    fimo_path: <Paht> a TXT file contains fimo output (*.fimo.txt).
    freq_thread: <str> minmum required binding frequence of TF.
    se_gene: <DataFrame> a DataFrame contains SEs regions and their     
        associated genes."""

    tf_bind = pd.read_csv(fimo_path, sep='\t').iloc[:, [0, 2]]

    tf_bind['binding_region'] = tf_bind['sequence_name'].map(
        lambda x: x.split(':')[0])

    tf_motif_df  = pd.read_table(TF_motif_dict_path, header=None)
    tf_motif_df.columns = ['motif_id', 'tf_name']

    tf_bind = pd.merge(tf_bind, tf_motif_df, on='motif_id', how='inner')

    se_gene = se_gene[['REGION_ID', 'name2']]

    tf_se_gene = pd.merge(
        tf_bind, se_gene,
        left_on='binding_region',
        right_on='REGION_ID',
        how='inner')

    tf_se_gene = tf_se_gene[~tf_se_gene.isna().any(axis=1)]

    freq = tf_se_gene.groupby(
        ['tf_name', 'name2'], as_index=False).count().iloc[:, :3]
    freq.columns = ['tf_name', 'name2', 'freq']

    freq = freq[freq['freq'] >= freq_thread]

    freq = pd.merge(
        freq, se_gene, 
        on='name2', how='inner'
    )

    freq = freq[['tf_name', 'REGION_ID', 'name2', 'freq']]

    return freq.copy()


def count_bind_freq_tss(fimo_path, refseq, odir, freq_thread=3, count_id_types='gn'):
    """Read in fimo output (*.fimo.txt), then count binding frequence of the TF.
    
    fimo_path: <Paht> a TXT file contains fimo output (*.fimo.txt).
    refseq: <DataFrame> ucsc gene annotation files.
    freq_thread: <str> minmum required binding frequence of TF.
    count_id_types: <str> 'gn' to count frequency by gene name,
        'tx' to count frequency by transcripts name.
    """
    id_colname = {
        'gn': 'gname',
        'tx': 'tx_id',
    }

    tf_tx_symbol = pd.read_csv(fimo_path, sep='\t').iloc[:, [0, 2]]

    tf_tx_symbol['tx_symbol'] = tf_tx_symbol['sequence_name'].map(
        lambda x: x.split(':')[0])
    tf_tx_symbol['bind_regions'] = tf_tx_symbol['sequence_name'].map(
        lambda x: x.split(':')[1].split('|')[0])

    tf_motif_df = pd.read_table(TF_motif_dict_path, header=None)
    tf_motif_df.columns = ['motif_id', 'tf_name']

    tf_tx_symbol = pd.merge(tf_tx_symbol, tf_motif_df, on='motif_id', how='inner')

    tx_symbol_tx = pd.DataFrame(
        list(
            chain(
                *pd.read_csv(
                    os.path.join(WorkDir, 'Tx_symbol_ids.dict.csv')
                ).drop_duplicates().apply(
                    lambda x: [[x[0], i] for i in x[1].split(',')], axis=1
                ).tolist()
            )
        ),
        columns = ['tx_symbol', 'tx_id']
    )

    tx_gn = refseq[['name', 'name2']]

    # annotate 'tx_symbol' by gene names
    tx_symbol_tx_gn = pd.merge(
        tx_symbol_tx, tx_gn,
        left_on='tx_id',
        right_on='name',
        how='inner'
    )
    tx_symbol_tx_gn.drop(['name'], axis=1)
    tx_symbol_tx_gn.rename({'name2': 'gname'}, axis=1, inplace=True)

    # relation of tf, tx_symbol, tx, genes
    tf_tx_symbol_tx_gn = pd.merge(
        tf_tx_symbol, tx_symbol_tx_gn,
        on='tx_symbol',
        how='inner'
    )

    tf_tx_symbol_tx_gn = tf_tx_symbol_tx_gn[
        ~tf_tx_symbol_tx_gn.isna().any(axis=1)]

    ofile = os.path.join(odir, 'tf_tss_txORg.raw.csv')
    tf_tx_symbol_tx_gn.to_csv(ofile, index=False)

    col_to_count = id_colname[count_id_types]

    freq = tf_tx_symbol_tx_gn.groupby(
        ['tf_name', col_to_count], as_index=False).count().iloc[:, :3]
    freq.columns = ['tf_name', col_to_count, 'freq']

    freq = freq[freq['freq'] >= freq_thread]

    freq = pd.merge(
        freq, tf_tx_symbol_tx_gn[['tx_symbol', 'bind_regions', col_to_count]],
        on=col_to_count, how='inner'
    )

    freq = freq[['tf_name', 'tx_symbol', 'bind_regions', col_to_count, 'freq']]

    return freq.drop_duplicates().copy()


def find_reliable_tf_se_gene(se_gene, se_fimo, odir, se_freq_thread=3):

    tf_se_gene_freq = count_bind_freq_se(se_fimo, se_gene, se_freq_thread)

    print('<TF SE Genes DataFrame:>\n', tf_se_gene_freq.head())

    ofile = os.path.join(odir, 'tf_se_gene_freq.csv')

    tf_se_gene_freq.to_csv(ofile, index=False)

    return tf_se_gene_freq.copy()


def find_reliable_tf_tss(fimo_path, refseq, odir, freq_thread=3, count_id_types='gn'):

    tf_tx_symbol_txORg_freq = count_bind_freq_tss(
        fimo_path, refseq, odir, freq_thread, count_id_types)

    print('<TF associated Genes or Tx freq:>\n',
          tf_tx_symbol_txORg_freq.head())

    ofile = os.path.join(odir, 'tf_tss_txORg_freq.csv')

    tf_tx_symbol_txORg_freq.to_csv(ofile, index=False)

    return tf_tx_symbol_txORg_freq.copy()


def find_reliable_se_gene(itf, se_gene, se_fimo, tss_fimo, odir, se_freq_thread=3, tss_freq_thread=3):

    se_reliable = count_bind_freq_se(se_fimo, se_freq_thread)
    se_reliable.columns = ['se_region', 'se_freq']
    se_reliable.to_csv(os.path.join(
        odir, 'SE_up_thread_regions.csv'), index=False)
    tss_reliable = count_bind_freq_se(tss_fimo, tss_freq_thread)
    tss_reliable.columns = ['g_region', 'g_freq']
    tss_reliable.to_csv(os.path.join(
        odir, 'TSS_up_thread_regions.csv'), index=False)

    se_gene = se_gene[
        (
            se_gene['se_name'].isin(se_reliable['se_region'])
        ) & (
            se_gene['g_name1'].isin(tss_reliable['g_region'])
            )
    ]

    se_gene['TF'] = itf

    se_gene = pd.merge(
        se_gene, se_reliable, left_on='se_name', right_on='se_region', how='outer')
    se_gene = pd.merge(
        se_gene, tss_reliable, left_on='g_name1', right_on='g_region', how='outer')

    se_gene = se_gene[~(se_gene.isna().any(axis=1))]

    #print(se_gene)

    se_gene = se_gene[
        ['TF', 'se_name', 'se_freq', 'g_name1', 'g_name2', 'g_freq']]

    ofile = os.path.join(odir, '{}_assignment.se.gene.csv'.format(itf))

    se_gene.to_csv(ofile, index=False)

    return ofile






