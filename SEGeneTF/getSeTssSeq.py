# -*- coding: utf-8 -*-
"""
Created on 2018-09-13 09:19:13
Last Modified on 2018-09-13 09:19:13

Get sequence from SE constitute and gene TSS region.

@Author: Ying Huang
"""
from collections import deque, namedtuple
import gzip
import os
import pandas as pd
from pybedtools.bedtool import BedTool

# self packages
from .checkFile import check_if_file_dup
from .regionTools import merge_overlaped_regions, get_tss


WorkDir = os.path.curdir
PWD = os.path.split(os.path.realpath(__file__))[0]
SplitChrDir = os.path.join(PWD, 'splitChr', 'gal5')
BED4_Format_Colname = ('chr', 'start', 'end', 'name')


def get_SE_constitute_regions(se_region, chip_sigal_region, extent=500):
    """Get SE constitute regions.
    
    se_region: <DataFrame> a dataframe contains SE regions in BED format.
    chip_sigal_region: <DataFrame> a dataframe contains ChIP-Seq signal regions
    in BED format.
    extent: <int> extent both site of each peaks."""

    extent = int(extent)

    se_bed = BedTool.from_dataframe(
        se_region.sort_values(['chr', 'start', 'end']))
    sig_bed = BedTool.from_dataframe(
        chip_sigal_region.sort_values(['chr', 'start', 'end']))

    # assign SEs to ChIP-Seq signals
    se_sig = se_bed.intersect(sig_bed, wa=True, wb=True).to_dataframe()

    # reset new name
    se_sig['name2'] = se_sig.iloc[:, [3,7]].apply(
        lambda x: '{}:{}'.format(str(x[0]), str(x[1])),
        axis=1
    )
    # format output
    se_sig = se_sig.iloc[:, [4,5,6,8]]
    se_sig.columns = BED4_Format_Colname

    # extent peaks
    se_sig['start'] = se_sig['start'].astype(int).map(
        lambda x: x - extent if x > extent else 0
    )
    se_sig['end'] = se_sig['end'].astype(int) + extent

    # merge overlaped extent peaks
    se_sig = merge_overlaped_regions(se_sig, 'simple_names_type')

    return se_sig.copy()


def get_tss_constitute_regions(annfile, chip_sigal_region, tss_up_extent=2.5e3, tss_dn_extent=2.5e3, extent=500):
    """Get ChIP-Seq signal regions in tss.
    
    annfile: <DataFrame> genes annotation file in a table separated format of   
        UCSC. Result of readFiles.read_annotation().
    chip_sigal_region: <DataFrame> a dataframe contains ChIP-Seq signal regions
        in BED format.
    tss_up_extent: <int> tss up stream extent length.
    tss_dn_extent: <int> tss down stream extent length.
    extent: <int> extent both site of each peaks."""

    tss_up_extent = int(tss_up_extent)
    tss_dn_extent = int(tss_dn_extent)
    extent = int(extent)

    # get TSS regions
    tss_df = get_tss(annfile, tss_up_extent, tss_dn_extent)

    tss_bed = BedTool.from_dataframe(
        tss_df.sort_values(['chr', 'start', 'end']))
    sig_bed = BedTool.from_dataframe(
        chip_sigal_region.sort_values(['chr', 'start', 'end']))

    # assign SEs to ChIP-Seq signals
    tss_sig = tss_bed.intersect(sig_bed, wa=True, wb=True).to_dataframe()

    # reset new name
    tss_sig['name2'] = tss_sig.iloc[:, [3, 7]].apply(
        lambda x: '{}:{}'.format(str(x[0]), str(x[1])),
        axis=1
    )
    # format output
    tss_sig = tss_sig.iloc[:, [4, 5, 6, 8]].drop_duplicates()
    tss_sig.columns = BED4_Format_Colname

    # extent peaks
    tss_sig['start'] = tss_sig['start'].astype(int).map(
        lambda x: x - extent if x > extent else 0
    )
    tss_sig['end'] = tss_sig['end'].astype(int) + extent

    # merge overlaped extent peaks
    tss_sig = merge_overlaped_regions(
        tss_sig, 'tss_names_type',
        os.path.join(WorkDir, 'Tx_symbol_ids.dict.csv')
    )

    return tss_sig.copy()


def read_fa(ipath):

    fasta = deque()

    with gzip.open(ipath, 'rt') as f:
        for line in f:
            if not line.startswith('>'):
                fasta.append(line.strip('\n'))

    return ''.join(fasta)


def search_seq(start, end, seq):
    """
    start: start site of seq, 0-base.
    end: end site of seq, 0-base.
    seq: DNA sequence of a chromsome.
    """

    assert isinstance(seq, str)
    assert isinstance(start, int) & isinstance(end, int)

    return seq[start: end + 1]


def get_seq(idir_fa, ibed, ofile):

    tmp_data = deque()

    #chroms = set(
    #    n[:-6] for n in os.listdir(idir_fa)
    #)
    # Get chromosome symbols intersection of BED file and FA files,
    # to make sure chromosome symbols are exist in both files.
    bed_chroms = set(ibed.chr.astype(str))
    fa_chroms = set([n.strip('.fa.gz') for n in os.listdir(idir_fa)])
    chroms = bed_chroms & fa_chroms
    unread_chroms = bed_chroms - fa_chroms
    if unread_chroms:
        print(
            '<Warnning:> some chromosomes are skiped:\n', unread_chroms,
            '\n\tintersect:\n', chroms,
            '\n\tchroms of bed\n:', bed_chroms,
            '\n\tchroms of fa:\n', fa_chroms,
        )


    for chrom in chroms:

        # Read separated FASTA file that generated by 'splitFA'.
        ifile_fa = os.path.join(idir_fa, str(chrom) + '.fa.gz')
        print("Reading FASTA file from '{}'...".format(ifile_fa))
        seq = read_fa(ifile_fa)

        for ch, start, end, name in ibed.loc[(ibed['chr'].astype(str) == chrom)].values:

            # Search sequence by BED file information.
            content = search_seq(start, end, seq)
            # Make header of searched seq to write out.
            header = ">{}|{}|{}|{}".format(
                name, ch, str(start), str(end))
            # Store searched sequence and its header in tmp variable
            tmp_data.append(header + '\n' + content + '\n')

            # Write out when tmp variable stored 1000 genes sequence
            # (header and content)
            if len(tmp_data) >= 1000:
                print("Output 1000 genes sequence...")
                with open(ofile, 'a') as f:
                    f.write(''.join(tmp_data))
                tmp_data.clear()

    print("Output final genes sequence.")
    if len(tmp_data) > 0:
            with open(ofile, 'a') as f:
                f.write(''.join(tmp_data))
            tmp_data.clear()


#==============#
#   main func  #
#==============#

def get_tss_fa(annfile, chip_sigal_region, odir=WorkDir, ifa_dir=SplitChrDir, tss_up_extent=2.5e3, tss_dn_extent=2.5e3, extent=500):
    """Extract H3K27ac ChIP-Seq regions in genes TSS regions."""
    
    tss_subpeaks = get_tss_constitute_regions(
        annfile, chip_sigal_region, tss_up_extent, tss_dn_extent, extent)

    tss_subpeaks_bed_ofile = os.path.join(odir, 'peaks_in_TSS.csv')
    tss_subpeaks.to_csv(tss_subpeaks_bed_ofile, index=False)

    tss_subpeaks_fa_ofile = os.path.join(odir, 'peaks_in_TSS.fa')
    # Make sure 'ofile' name hasn't been used before.
    tss_subpeaks_fa_ofile = check_if_file_dup(
        tss_subpeaks_fa_ofile, if_rmfile=True)

    print('<Get seq of peaks in TSS {}-{} ...>'.format(tss_up_extent, tss_dn_extent))
    get_seq(ifa_dir, tss_subpeaks, tss_subpeaks_fa_ofile)

    Ofiles = namedtuple('Ofiles', ['TSS_BED', 'TSS_FA'])

    return Ofiles(tss_subpeaks_bed_ofile, tss_subpeaks_fa_ofile)


def get_se_fa(se_region, chip_sigal_region, odir=WorkDir, ifa_dir=SplitChrDir, extent=500):

    se_subpeaks = get_SE_constitute_regions(se_region, chip_sigal_region, extent)

    se_subpeaks_bed_ofile = os.path.join(odir, 'peaks_in_SE.csv')
    se_subpeaks.to_csv(se_subpeaks_bed_ofile, index=False)

    se_subpeaks_fa_ofile = os.path.join(odir, 'peaks_in_SE.fa')
    # Make sure 'ofile' name hasn't been used before.
    se_subpeaks_fa_ofile = check_if_file_dup(se_subpeaks_fa_ofile, if_rmfile=True)

    print('<Get seq of peaks in SE ...>')
    get_seq(ifa_dir, se_subpeaks, se_subpeaks_fa_ofile)

    Ofiles = namedtuple('Ofiles', ['SE_BED', 'SE_FA'])

    return Ofiles(se_subpeaks_bed_ofile, se_subpeaks_fa_ofile)



