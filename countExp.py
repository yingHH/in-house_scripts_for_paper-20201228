# -*- coding: utf-8 -*-
"""
Author: Ying Huang

Date: 2020-04-08 00:27:36
Last Modified by: Ying Huang
Last Modified time: 2020-04-08 00:27:36

Description: 
Calculate gene expression
"""
import pandas as pd
from collections import deque
import re
from pybedtools.bedtool import BedTool


def ncbi_exon_lst(igff):
    """Make list of gene refseq ids and Gene names, Gene ids.
    
    igff: <Path> path to NCBI RefSeq file in GFF foramt."""

    exon_list = deque()
    
    with open(igff, 'rt') as f:
        for line in f:
            fields = line.split('\t')
            if (not line.startswith('#')) and (fields[2] == 'exon'):
                exon = re.search(r'ID=([^;]+);', line)
                gnames = re.search(r'gene=([^;]+);', line)
                start = fields[3]
                end = fields[4]
                chrom = fields[0]
                if exon and gnames:
                    exon_list.append([chrom, start, end, exon.group(1), gnames.group(1)])
    # make DataFrame, columns are 'mask_id', 'gene_name', 'general_id'
    exon_id_df = pd.DataFrame(list(exon_list))
    exon_id_df.columns = ['chr', 'start', 'end', 'exons', 'genes']

    return exon_id_df.copy()

def ens_exon_lst(igff):
    """Make list of gene ENS ids and Gene names, Gene ids.
    
    igff: <Path> path to ENSEMBL GFF3 file in GFF foramt."""

    exon_list = deque()
    
    with open(igff, 'rt') as f:
        for line in f:
            fields = line.split('\t')
            if (not line.startswith('#')) and (fields[2] == 'exon'):
                exon = re.search(r'ID=([^;]+);', line)
                gnames = re.search(r'gene=([^;]+);', line)
                start = fields[3]
                end = fields[4]
                chrom = fields[0]
                if exon and gnames:
                    exon_list.append([chrom, start, end, exon.group(1), gnames.group(1)])
    # make DataFrame, columns are 'mask_id', 'gene_name', 'general_id'
    exon_id_df = pd.DataFrame(list(exon_list))
    exon_id_df.columns = ['chr', 'start', 'end', 'exons', 'genes']

    return exon_id_df.copy()


def gene_exons_len(exon_df):
    
    df = exon_df[['genes','start', 'end']]
    df.iloc[:, 1:] = df.iloc[:, 1:].astype(int)
    df.sort_values(df.columns.tolist(), inplace=True)
    print('<Merge exons regions by bedtools> ...')
    exon = BedTool.from_dataframe(df)
    exons_pos = exon.merge().to_dataframe()
    
    exons_pos['exons_len'] = exons_pos['end'] - exons_pos['start'] + 1
    print('<Calculate gene exons length> ...')
    res = exons_pos[['chrom', 'exons_len']].groupby('chrom').agg('sum')
     
    return res.copy()


def cal_exonlen(igff):
    print('<List exons of>: ', igff)
    exon_lst = ncbi_exon_lst(igff)
    gexon_len = gene_exons_len(exon_lst)
    gexon_len.index.name = 'genes'
    return gexon_len.copy()


def rpkm(df, igff):
    """
    Calculate RPKM.

    df: <DataFrame> first column is gene names, other column names are sample names, contents are reads count.
    igff: <Path> GFF3 file from Refseq NCBI.
    """
    gexon_len = cal_exonlen(igff)
    df.drop_duplicates(df.columns[0], inplace=True)
    df.set_index(df.columns[0], drop=True, inplace=True)
    
    print('''
<Gene sets information>:
    Number of 'GFF3' genes: {}
    All detected genes: {}
    Detected genes in 'hg38 refseq': {}
    '''.format(
            len(gexon_len), len(df), len(set(gexon_len.index.tolist()) & set(df.index.tolist()))
        )
    )

    detected_genes = list(set(gexon_len.index.tolist()) & set(df.index.tolist()))

    df = df.loc[detected_genes, :]
    print('<Shape of DataFrame to calculated:> ', df.shape)
    #count total reads
    total_reads = df.sum(axis=0)
    # trans total reads to millions (10^-6)
    total_reads = total_reads.map(lambda x: float(x) * 10 ** -6)

    # get exon len
    gexon_len = gexon_len.loc[detected_genes, 'exons_len']
    # trans exon len to KB (10^-3)
    gexon_len = gexon_len.map(lambda x: float(x) * 10 ** -3)

    print('<Calculate RPKM>')
    rpkm = df.astype(float) / total_reads / gexon_len.values.reshape(len(gexon_len), 1)
    print('done')

    return rpkm

    