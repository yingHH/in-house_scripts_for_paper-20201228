# -*- coding: utf-8 -*-
"""
Created on 2018-11-30 23:19:59
Last Modified on 2018-11-30 23:19:59

Run fimo to find TF binding site to peaks in SE or TSS.

@Author: Ying Huang
"""
from collections import namedtuple, defaultdict, Iterable
import os
import pandas as pd


PWD = os.path.split(os.path.realpath(__file__))[0]
Motif_db_path = os.path.join(PWD, 'motif_db', 'VertebratePWMs.txt')
TF_motif_dict_path = os.path.join(PWD, 'motif_db', 'MotifDictionary.txt')

def mk_tf_motif_dict(ipath):

    tf_motif_dict = defaultdict(list)
    tf_motif_df = pd.read_csv(
        ipath, header=None, sep='\t').iloc[:, [1,0]]
    
    tf_motif_df[1] = tf_motif_df[1].map(lambda x: x.upper())

    for k, v in tf_motif_df.values:
        tf_motif_dict[k].append(v)

    return tf_motif_dict


def fimo_bg(ifa, odir, oname):
    """Make fimo backgroud sequence by sequence in SE or TSS.
    
    ifa: <Path> FASTA file of peaks in SE or TSS.
    odir: <Path> ouput directory.
    oname: <str> ouput name."""

    oname = '{}_bg.meme'.format(oname)
    ofile = os.path.join(odir, oname)

    cmd = \
    "/usr/local/bin/meme/libexec/meme-5.0.2/fasta-get-markov \
    -m 1 \
    {} \
    {}".format(ifa, ofile)

    os.system(cmd)

    return ofile


def run_fimo(itfs, ifa, ibg, odir, oname, TFmotifDict):
    """Find TF binding regions by sequence in SE or TSS.
    
    itfs: <iter> tf names.
    ifa: <Path> FASTA file of peaks in SE or TSS.
    odir: <Path> ouput directory.
    oname: <str> ouput name."""

    ofile = os.path.join(odir, '{}.fimo.txt'.format(oname))

    #print('<itfs:>', itfs)

    print('<Input TFs:>', [itf for itf in itfs])

    motif_cmd = ' '.join([
        "--motif '{}'".format(motif) for itf in itfs \
        for motif in TFmotifDict[str(itf).upper()]
    ])

    cmd = \
    "/usr/local/bin/MEME/fimo \
    {} \
    -verbosity 1 \
    -text \
    -oc {} \
    --bgfile {} \
    {} \
    {} \
    > {}".format(
        motif_cmd,
        odir,
        ibg,
        Motif_db_path,
        ifa,
        ofile
    )

    res = os.system(cmd)

    print('<Fimo cmd:>\n', cmd, '\n', res)

    return ofile


#
#   main func
#

def fimo(itfs, ifa, odir, oname):

    itfs = [itfs] if not isinstance(itfs, Iterable) else itfs

    TFmotifDict = mk_tf_motif_dict(TF_motif_dict_path)

    Ofiles = namedtuple('Ofiles', ['fimo', 'bg'])

    bg = fimo_bg(ifa, odir, oname)
    fimo = run_fimo(itfs, ifa, bg, odir, oname, TFmotifDict)

    return Ofiles(fimo, bg)
