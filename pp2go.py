#!/usr/bin/env python3

from typing import Dict, List
import pandas as pd
import os
import logging

import aiohttp
import asyncio
import pdb
import sys

from phylogenetic_profiling import PhylogeneticProfiling

__author__ = "Geovanny Risco"
__email__ = "geovanny.risco@bsc.es"
__version__ = "0.1"



__DESCRIPTION__ = ""


def intersectSwissProt(uniprotIDS: set, swiss_prot_ids: os.path):
    with open(swiss_prot_ids, "r") as input:
        sp_ids = input.read().split()
    sp_ids = list(map(lambda line: line.split(","), sp_ids))
    sp_ids_exploded = [item for packed_elements in sp_ids for item in packed_elements]
    return uniprotIDS.intersection(set(sp_ids_exploded))

def filterBySwissProt(df: pd.DataFrame, onColumns: List[str], swiss_prot_ids: os.path = "drive/MyDrive/TFG/swiss_prot_ids.csv"):
    for column in onColumns:
        swiss_prot_proteins = intersectSwissProt(set(df[column].unique()), swiss_prot_ids)
        df = df[df[column].isin(swiss_prot_proteins)]
    return df

def main():
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] -- %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stdout
    )    

    mtp_orthologs = "drive/MyDrive/TFG/QfO_input.tsv"

    all_orthologs = pd.read_table(mtp_orthologs, names=["uniprotid1", "uniprotid2"], dtype="string")
    orthologs_swiss_prot_only = filterBySwissProt(all_orthologs, onColumns=all_orthologs.columns)
    pp = PhylogeneticProfiling(orthologs_swiss_prot_only, onSpecies=["9606"])
    pp_matrix = pp.computeCountsMatrix()
    pp_matrix.to_csv("drive/MyDrive/TFG/orthologs_counts_matrix_v2.tsv", sep="\t", index=True, header=True)
    
    human_orthologs = ""
    genes2GOterms = None
    # TODO: train_model()
    pass

if __name__ == "__main__":
    main()
      
