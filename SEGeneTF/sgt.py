# -*- coding: utf-8 -*-
"""
Created on 2018-11-30 17:26:21
Last Modified on 2018-11-30 17:26:21

Main function to find TFs associated super enhancers
and genes regulated by TFs.

@Author: Ying Huang
"""
from collections import namedtuple, Iterable
import os
import pandas as pd

# self packages
from .readFiles import *
from .seGene import assign_se_to_genes
from .getSeTssSeq import get_tss_fa, get_se_fa
from .runfimo import fimo
from .TFSEGene import find_reliable_se_gene, find_reliable_tf_se_gene, find_reliable_tf_tss


WorkDir = os.path.curdir
PWD = os.path.split(os.path.realpath(__file__))[0]
SplitChrDir = os.path.join(PWD, 'splitChr', 'gal5')
BED4_Format_Colname = ('chr', 'start', 'end', 'name')
BED6_Format_Colname = ('chr', 'start', 'end', 'name', 'score', 'strand')

SP_Content = namedtuple('SP_Content', ['TF_db'])
SP_DICT = {
    'gal5': SP_Content(os.path.join(PWD, 'TF_db', 'TFlist_NMid_gal5.txt')),
}


def tf_tss_rna(ann_path, chip_sig_path, tfs_input=None, crc_tfs_path=None, odir=WorkDir, ifa_dir=SplitChrDir, if_rm_chr=True, exp_genes=None, genes_exp_path=None, exp_id_type='tx', tss_up_extent=2.5e3, tss_dn_extent=2.5e3, extent=500, tss_freq_thread=3, count_id_types='gn'):

    Res = namedtuple(
        'Res',
        [
            'Refseq', 'ExpGene', 'ChipSig',
            'tss_subpeaks_bed_ofile', 'tss_subpeaks_fa_ofile', 
            'tss_fimo', 'tss_bg',
            'tf_tx_symbol_txORg_freq',
        ]
    )

    if not os.path.isdir(odir):
        os.mkdir(odir)

    exp_csv = exp_lst = set()
    if genes_exp_path is not None:
        tmp_exp = pd.read_csv(genes_exp_path, header=None)[0]
        tmp_exp.to_csv(
            os.path.join(odir, 'expressed_genes.csv'),
            index=False, header=None)
        exp_csv = set(tmp_exp)

    if exp_genes is not None:
        assert (not isinstance(exp_genes, str)
                ) and isinstance(exp_genes, Iterable)
        exp_lst = set(exp_genes)

    exp = tuple(exp_csv | exp_lst)
    exp = exp if bool(exp) else None

    # read input
    refseq = read_annotation(ann_path, if_rm_chr, exp, exp_id_type)

    #tf_db = read_tf_db(SP_DICT[species].TF_db)
    tfs = crc_tfs = set()
    if tfs_input:
        tfs = read_tfs(tfs_input)
    if crc_tfs_path and os.path.exists(crc_tfs_path):
        crc_tfs = read_crc_tf(crc_tfs_path)
    if tfs_input is None and crc_tfs_path is None:
        raise Exception(
            '<ERR:> at least one value should be passed to "tfs_input" or "crc_tfs_path".'
        )
    cand_tf = tfs | crc_tfs
    print('<Candidate TFs:> ', cand_tf)

    chip_sig = read_chip_sigal_regions(chip_sig_path)

    # get sequence of Peaks in TSS region
    tss_subpeaks_bed_ofile, tss_subpeaks_fa_ofile = get_tss_fa(
        refseq, chip_sig, odir, ifa_dir, tss_up_extent, tss_dn_extent, extent)

    # run fimo
    tss_fimo, tss_bg = fimo(cand_tf, tss_subpeaks_fa_ofile, odir, 'tf_tss')

    # get reliable Se gene assignment
    ofile = find_reliable_tf_tss(
        tss_fimo, refseq, odir, tss_freq_thread, count_id_types)

    res = Res(
        refseq, exp, chip_sig,
        tss_subpeaks_bed_ofile, tss_subpeaks_fa_ofile,
        tss_fimo, tss_bg,
        ofile,
    )

    return res


def se_gene_tf(se_path, ann_path, chip_sig_path, tfs_input=None, crc_tfs_path=None, odir=WorkDir, ifa_dir=SplitChrDir, if_rm_chr=True,  exp_genes=None, genes_exp_path=None, extent=500, se_freq_thread=3):

    Res = namedtuple(
        'Res',
        [
            'SE', 'Refseq', 'ExpGene', 'ChipSig', 'SEGene',
            'se_bed', 'gene_bed', 'se_subpeaks_bed_ofile', 'se_subpeaks_fa_ofile', 'se_fimo', 'se_bg',
            'TF_SE_Genes',
        ]
    )

    if not os.path.isdir(odir):
        os.mkdir(odir)

    # read input
    se = read_se(se_path)
    refseq = read_annotation(ann_path, if_rm_chr)

    exp_csv = exp_lst = set()
    if genes_exp_path is not None:
        tmp_exp = pd.read_csv(genes_exp_path, header=None)[0]
        tmp_exp.to_csv(
            os.path.join(odir, 'expressed_genes.csv'), 
            index=False, header=None)
        exp_csv = set(tmp_exp)

    if exp_genes is not None:
        assert (not isinstance(exp_genes, str)) and isinstance(exp_genes, Iterable)
        exp_lst = set(exp_genes)

    exp = tuple(exp_csv | exp_lst)

    #tf_db = read_tf_db(SP_DICT[species].TF_db)
    tfs = crc_tfs = set()
    if tfs_input:
        tfs = read_tfs(tfs_input)
    if crc_tfs_path and os.path.exists(crc_tfs_path):
        crc_tfs= read_crc_tf(crc_tfs_path)
    if tfs_input is None and crc_tfs_path is None:
        raise Exception(
            '<ERR:> at least one value should be passed to "tfs_input" or "crc_tfs_path".'
        )
    cand_tf = tfs | crc_tfs
    print('<Candidate TFs:> ', cand_tf)

    chip_sig = read_chip_sigal_regions(chip_sig_path)

    # assign SEs to genes
    se_genes = assign_se_to_genes(se, refseq, exp)

    # make SE region in BED format from SE_Gene result
    se_bed = se_genes.iloc[:, :4]
    se_bed.columns = BED4_Format_Colname
    print('<SE BED:>\n', se_bed.head())

    # make Gene region in BED format from SE_Gene result
    gene_bed = se_genes.iloc[:, 4:]
    #gene_bed.columns = BED4_Format_Colname
    print('<Gene BED:>\n', gene_bed.head())

    # get sequence of SE and TSS region
    se_subpeaks_bed_ofile, se_subpeaks_fa_ofile = get_se_fa(
        se_bed, chip_sig, odir, ifa_dir, extent)

    # run fimo
    se_fimo, se_bg = fimo(cand_tf, se_subpeaks_fa_ofile, odir, 'tf_se')

    # get reliable Se gene assignment
    ofile = find_reliable_tf_se_gene(
        se_genes, se_fimo, odir, se_freq_thread,)

    res = Res(
        se, refseq, exp, chip_sig, se_genes,
        se_bed, gene_bed, se_subpeaks_bed_ofile, 
        se_subpeaks_fa_ofile, se_fimo, se_bg,
        ofile,
    )

    return res
