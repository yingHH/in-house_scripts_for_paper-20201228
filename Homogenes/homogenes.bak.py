# -*- coding: utf-8 -*-
"""
Created on 2018-11-29 15:48:26
Last Modified on 2018-11-29 15:48:26

The main script to make homogenes list.

@Author: Ying Huang
"""
import argparse
import os
import re
import pandas as pd
import datatable as dt
from collections import deque, namedtuple
from concurrent import futures
import multiprocessing
from functools import reduce

# self library
from Homogenes.find_homogenes import find_homogenes
from Homogenes.add_refseq_ids import add_gene_ids
from Homogenes.speices_name_id_lst import SpNameIdLst


PWD = os.path.split(os.path.realpath(__file__))[0]
HOMO_DB = os.path.join(PWD, 'homo_db', 'biuld68')
Species_DB_ID = pd.read_csv(os.path.join(HOMO_DB ,'taxid_taxname'), sep='\t', header=None)


def convert_sp_idx_to_id(sp_name_id_obj, each_sp_idx):
    splst = sp_name_id_obj

    if each_sp_idx in splst._taxid_taxname_df['taxonomy_name'].tolist():
            return splst.search(each_sp_idx)
    else:
        if each_sp_idx in splst._taxid_taxname_df['taxonomy_id'].astype(str).tolist():
            return int(each_sp_idx)
        else:
            raise Exception("<ERR:> species symbol '{}' not found.".format(each_sp_idx))


def merge_same_gnames_in_refseq(gffs):
    """
    Gene names that are same in NCBI RefSeq across species are regarded as homologous genes.

    Required:
    gffs: <list> a list of each species gff3 PATH of NCBI.
    """
    spec_gff_lst = []
    def add_gff_to_lst(gff):
        gnames_deque = deque()
        with open(gff, 'rt') as f:
            for line in f:
                if re.search(r'ID=gene', line):
                    names = re.search(r'Name=([^;]+);', line)
                    gnames_deque.append(names.group(1))
        gnames_lst = list(set(gnames_deque))
        gnames_dict = {i.upper(): i for i in gnames_lst}
        df = dt.Frame([gnames_dict])
        df[0, 'SP_GFF'] = gff
        spec_gff_lst.append(df)

        return 'read in genes: {}'.format(str(len(gnames_lst)))

    isCancel = False
    with futures.ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        toDo = []
        for order, gff in enumerate(gffs):
            future = executor.submit(add_gff_to_lst, gff)
            toDo.append(future)
            msg = '<Reading in {}>: \n\t{}'
            print(msg.format(gff, future))
        try:
            results = []
            for future in futures.as_completed(toDo):
                res = future.result()
                msg = '<{} result>: \n\t{!r}'
                print(msg.format(future, res))
                results.append(res)
        except KeyboardInterrupt:
            isCancel = True
            for future in futures:
                future.cancel()

    print('<Merge all gene names> ...')
    # get common genes across species
    comm_gnames = list(reduce(set.intersection, [set(df.names) for df in spec_gff_lst]))
    # remove uncommon genes of each Dataframe
    spec_gff_lst_fillter = [df[:, comm_gnames] for df in spec_gff_lst]
    # merge all DataFrame
    res = dt.rbind(spec_gff_lst_fillter, force=True).to_pandas()

    return res.copy()

    
#==========#
#   Args   #
#==========#

def cmd():

    msg = """
    make homogenes list by NCBI homogenes files.
    (download from: ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current)
    """

    parser = argparse.ArgumentParser(description=msg)

    # Required arguements

    parser.add_argument(
        'sp_idx',
        nargs='+',
        action='store',
        help="<list> species IDs or Names.(eg. 10090 'Gallus gallus' 'Mus musculus' 9031). Note: species Names should surround by quotation mark. Species ids and names can be listed by command '-lst'."
    )

    # Optional arguements

    parser.add_argument(
        '-of', '--ofile',
        dest='ofile',
        type=str,
        default='homogens_out.csv',
        action='store',
        required=False,
        help="<Path> ofile path."
    )

    parser.add_argument(
        '-db', '--homo_db',
        dest='homo_db',
        type=str,
        default=HOMO_DB,
        action='store',
        required=False,
        help="<Path> database directory contained files 'homologene.data, taxid_taxname, all_proteins.data'. All of the file can be download    from 'ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current'. Default is {}, which can be provide by cutomers.".format(HOMO_DB)
    )

    args = parser.parse_args()

    return args


#==========#
#   main   #
#==========#

def homogenes(sp_idx, gffs, homo_db=HOMO_DB):
    """
    Get homologous genes from NCBI Homologous gene databese ('ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current').

    Required:
    sp_idx: <list> species IDs or Names.(eg. 10090 'Gallus gallus' 'Mus musculus' 9031). 
        Note: species Names should surround by quotation mark.
    gffs: <list> a list of each species gff3 file of NCBI. Note: gffs order must be consistent with sp_idx order.
    homo_db: <Path> database directory contained files 'homologene.data, taxid_taxname, all_proteins.data'. 
        All of the file can be download from 'ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current'.
    """
    print('<!!! Note: gffs order must be consistent with sp_idx order !!!> ')
    homo_db_path = homo_db
    taxid_taxname_path = os.path.join(homo_db_path, 'taxid_taxname')

    splst = SpNameIdLst(taxid_taxname_path)

    sp_ids = [convert_sp_idx_to_id(splst, i) for i in sp_idx]

    print('<Finding homogenes of:> {}'.format(','.join(sp_idx)))
    homogenes_df = find_homogenes(homo_db_path, sp_ids)
    print('<Adding refseq gene ids ...>')
    homogenes_df = add_gene_ids(homo_db_path, homogenes_df, sp_ids)

    # convert species ids to names
    new_col = [
        '_'.join(
            n.split('_')[:-1] + [splst.search(int(n.split('_')[-1]))]
        ) for n in homogenes_df.columns
    ]

    homogenes_df.columns = new_col

    return homogenes_df

def main():
    args = cmd()
    res = homogenes(args.sp_idx, args.homo_db)

    ofile = args.ofile
    n = 0
    while os.path.exists(ofile):
        ofile = '{}{}'.format(ofile, n)
        n+=1

    res.to_csv(ofile, index=False)
    print('<Write to file:> {}'.format(ofile))
    

if __name__ == '__main__':
    print(
        """species IDs and Names:
        10090	Mus musculus
        10116	Rattus norvegicus
        28985	Kluyveromyces lactis
        318829	Magnaporthe oryzae
        33169	Eremothecium gossypii
        3702	Arabidopsis thaliana
        4530	Oryza sativa
        4896	Schizosaccharomyces pombe
        4932	Saccharomyces cerevisiae
        5141	Neurospora crassa
        6239	Caenorhabditis elegans
        7165	Anopheles gambiae
        7227	Drosophila melanogaster
        7955	Danio rerio
        8364	Xenopus (Silurana) tropicalis
        9031	Gallus gallus
        9544	Macaca mulatta
        9598	Pan troglodytes
        9606	Homo sapiens
        9615	Canis lupus familiaris
        9913	Bos taurus"""
    )
    main()

