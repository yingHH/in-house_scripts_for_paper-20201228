# -*- coding: utf-8 -*-
"""
Created on 2018-11-30 13:07:36
Last Modified on 2018-11-30 13:07:36

Assign SEs to closest genes.

@Author: Ying Huang
"""
import pandas as pd
from pybedtools.bedtool import BedTool


def assign_se_to_genes(se, refseq, exp_genes=None, up_extend=1e6, dn_extend=1e6):
    """Assign SEs to closest genes within 1Mb around SEs.
    
    se: <DataFrame> result of function readFiles.read_se().
    refseq: <DataFrame> result of function readFiles.read_annotation().
    exp_genes: <tuple> a list of expressed genes.
    up_extend: <int> length SEs upstream extended region.
    dn_extend: <int> length SEs downstream extended region.
    """

    # only use exp_genesressed genes for analysis below
    #print('<Refseq len:> ', len(refseq))
    if exp_genes is not None:
        refseq = refseq[refseq['name'].isin(exp_genes)]
    #print('<Refseq len:> ', len(refseq))

    #print(se.head())
    print(refseq.head())

    se_midpoint = pd.DataFrame()
    tss = pd.DataFrame()

    # make SE midpoint DataFrame
    se_midpoint['CHROM'] = se['CHROM']
    se_midpoint['START'] = se[['START', 'STOP']].astype(int).apply(
        lambda x: int((x['START'] + x['STOP']) / 2),
        axis=1,
    )
    se_midpoint['STOP'] = se[['START', 'STOP']].astype(int).apply(
        lambda x: int((x['START'] + x['STOP']) / 2),
        axis=1,
    )
    se_midpoint['REGION_ID'] = se['REGION_ID']
    print(se_midpoint.head())

    # record extented SEs boundaries
    se_midpoint['up_boundary'] = se['START'].astype(int) - up_extend
    se_midpoint['dn_boundary'] = se['STOP'].astype(int) + dn_extend

    # make refseq TSS DataFrame
    tss['chrom'] = refseq['chrom']
    # Get the TSS position, use the start site of annotation file 
    # if strand id '+', else use the end site.
    tss_pos = refseq.apply(
        lambda x: x[3] if x[2] == '+' else x[4],
        axis=1,    
    )

    tss['tssStart'] = tss_pos
    tss['tssEnd'] = tss_pos
    tss['name'] = refseq['name']
    print(tss.head())

    se_midpoint_bed = BedTool.from_dataframe(
        se_midpoint.sort_values(
            ['CHROM', 'START', 'STOP'])
    )
    tss_bed = BedTool.from_dataframe(
        tss.sort_values(['chrom', 'tssStart', 'tssEnd'])
    )

    close_res = se_midpoint_bed.closest(tss_bed, d=True).to_dataframe()

    # only save closest TSS in extented SEs regions
    print('<Raw SE-genes:>\n', len(close_res))
    print(close_res.head())
    mask = close_res.apply(
        lambda x: bool(int(x[4]) < int(x[1]) + int(x[10]) < int(x[5])),
        axis=1,
    )
    close_res = close_res[mask]

    print('<closest TSS in extend SEs:>\n', len(close_res))
    print(
        '<Duplicated SE:>\n', 
        close_res[close_res.duplicated([close_res.columns[3]])]
    )

    # remove unassigned rows
    mask = (close_res == -1)
    close_res = close_res[~(mask.any(axis=1))]

    # get SE ids and refseq gene ids
    close_rg_gene_ids = close_res.iloc[:, [3, 9]]
    close_rg_gene_ids.columns = ['SE_ids', 'refseq_ids']
    print(close_rg_gene_ids.head())

    ## add region information of SE and genes
    # add SEs information
    close_rg_gene_ids = pd.merge(
        close_rg_gene_ids, se,
        left_on='SE_ids', right_on='REGION_ID',
        how='outer',
    )
    # remove unmatched rows
    close_rg_gene_ids = close_rg_gene_ids[
        ~(close_rg_gene_ids.isna().any(axis=1))
    ]
    # add genes information
    print('<Refseq len:> ', len(refseq))
    close_rg_gene_ids = pd.merge(
        close_rg_gene_ids, refseq,
        left_on='refseq_ids', right_on='name',
        how='outer',
    )
    # remove unmatched rows
    close_rg_gene_ids = close_rg_gene_ids[
        ~(close_rg_gene_ids.isna().any(axis=1))
    ]

    # remove first two columns
    close_rg_gene_ids.drop(['SE_ids', 'refseq_ids'], inplace=True, axis=1)

    # Format output
    close_rg_gene_ids[['START', 'STOP', 'txStart', 'txEnd']] = \
    close_rg_gene_ids[['START', 'STOP', 'txStart', 'txEnd']].astype(int)

    close_rg_gene_ids[
        ['CHROM', 'REGION_ID', 'chrom', 'name', 'name2', 'strand']
    ] = close_rg_gene_ids[
        ['CHROM', 'REGION_ID', 'chrom', 'name', 'name2', 'strand']
    ].astype(str)

    close_rg_gene_ids

    close_rg_gene_ids = close_rg_gene_ids[
        ['CHROM', 'START', 'STOP', 'REGION_ID',
        'chrom', 'txStart', 'txEnd', 'name', 'name2', 'strand']
    ]
    
    return close_rg_gene_ids


def assign_rg_to_genes(rg, refseq, exp=pd.DataFrame(), up_extend=1e6, dn_extend=1e6):
    """Assign RGs to closest genes within 1Mb around RGs.
    
    rg: <DataFrame> a region contained "'CHROM', 'START', 'STOP', 'REGION_ID'" 
        four columns.
    refseq: <DataFrame> result of function readFiles.read_annotation().
    exp: <DataFrame> result of function readFiles.read_exp_genes(), if empty,
        no expression genes will be used.
    up_extend: <int> length RGs upstream extended region.
    dn_extend: <int> length RGs downstream extended region.
    """

    rg.columns = ['CHROM', 'START', 'STOP', 'REGION_ID']
    #print(rg.head())

    # only use expressed genes for analysis below
    #print('<Refseq len:> ', len(refseq))
    if not exp.empty:
        refseq = refseq[refseq['name'].isin(exp['ids'])]
    #print('<Refseq len:> ', len(refseq))

    #print(rg.head())
    print('<RefSeq:>\n', refseq.head())

    #rg_midpoint = pd.DataFrame()
    #tss = pd.DataFrame()

    # make RG midpoint DataFrame
    rg_midpoint = rg[['CHROM']].copy()
    rg_midpoint['START'] = rg[['START', 'STOP']].astype(int).apply(
        lambda x: int(x.sum() / 2),
        axis=1,
    )
    rg_midpoint['STOP'] = rg[['START', 'STOP']].astype(int).apply(
        lambda x: int(x.sum() / 2),
        axis=1,
    )
    rg_midpoint['REGION_ID'] = rg['REGION_ID']
    print('<rg_midpoint:>\n', rg_midpoint.head())

    # record extented RGs boundaries
    rg_midpoint['up_boundary'] = rg['START'].astype(int) - up_extend
    rg_midpoint['dn_boundary'] = rg['STOP'].astype(int) + dn_extend

    # make refseq TSS DataFrame
    tss = refseq[['chrom']].copy()
    # Get the TSS position, use the start site of annotation file
    # if strand id '+', else use the end site.
    tss_pos = refseq.apply(
        lambda x: x[3] if x[2] == '+' else x[4],
        axis=1,
    )

    tss['tssStart'] = tss_pos
    tss['tssEnd'] = tss_pos
    tss['name'] = refseq['name']
    print('<tss:>\n', tss.head())

    #print('<rg_midpoint:>\n', rg_midpoint.head())

    rg_midpoint_bed = BedTool.from_dataframe(
        rg_midpoint.sort_values(
            ['CHROM', 'START', 'STOP'])
    )
    tss_bed = BedTool.from_dataframe(
        tss.sort_values(['chrom', 'tssStart', 'tssEnd'])
    )

    close_res = rg_midpoint_bed.closest(tss_bed, d=True).to_dataframe()

    # only save closest TSS in extented RGs regions
    print('<Raw RG-genes:>\n', len(close_res))
    print(close_res.head())
    mask = close_res.apply(
        lambda x: bool(int(x[4]) < int(x[1]) + int(x[10]) < int(x[5])),
        axis=1,
    )
    close_res = close_res[mask]

    print('<closest TSS in extend RGs:>\n', len(close_res))
    print(
        '<Duplicated RG:>\n',
        close_res[close_res.duplicated([close_res.columns[3]])]
    )

    # remove unassigned rows
    mask = (close_res == -1)
    close_res = close_res[~(mask.any(axis=1))]

    # get RG ids and refseq gene ids
    close_rg_gene_ids = close_res.iloc[:, [3, 9]]
    close_rg_gene_ids.columns = ['RG_ids', 'refseq_ids']
    print(close_rg_gene_ids.head())

    ## add region information of RG and genes
    # add RGs information
    close_rg_gene_ids = pd.merge(
        close_rg_gene_ids, rg,
        left_on='RG_ids', right_on='REGION_ID',
        how='outer',
    )
    # remove unmatched rows
    close_rg_gene_ids = close_rg_gene_ids[
        ~(close_rg_gene_ids.isna().any(axis=1))
    ]
    # add genes information
    print('<Refseq len:> ', len(refseq))
    close_rg_gene_ids = pd.merge(
        close_rg_gene_ids, refseq,
        left_on='refseq_ids', right_on='name',
        how='outer',
    )
    # remove unmatched rows
    close_rg_gene_ids = close_rg_gene_ids[
        ~(close_rg_gene_ids.isna().any(axis=1))
    ]

    # remove first two columns
    close_rg_gene_ids.drop(['RG_ids', 'refseq_ids'], inplace=True, axis=1)

    # Format output
    close_rg_gene_ids[['START', 'STOP', 'txStart', 'txEnd']] = \
        close_rg_gene_ids[['START', 'STOP', 'txStart', 'txEnd']].astype(int)

    close_rg_gene_ids[
        ['CHROM', 'REGION_ID', 'chrom', 'name', 'name2', 'strand']
    ] = close_rg_gene_ids[
        ['CHROM', 'REGION_ID', 'chrom', 'name', 'name2', 'strand']
    ].astype(str)

    close_rg_gene_ids

    close_rg_gene_ids = close_rg_gene_ids[
        ['CHROM', 'START', 'STOP', 'REGION_ID',
         'chrom', 'txStart', 'txEnd', 'name', 'name2', 'strand']
    ]

    close_rg_gene_ids.drop_duplicates(
        ['REGION_ID', 'name2'],
        inplace=True,
    )

    return close_rg_gene_ids
