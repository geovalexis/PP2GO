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
from gene_ontology import GeneOntology, GO_Aspects, GO_EvidenceCodes

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

    # Define orthologs dataset
    all_orthologs = pd.read_table(mtp_orthologs, names=["uniprotid1", "uniprotid2"], dtype="string")
    orthologs_swiss_prot_only = filterBySwissProt(all_orthologs, onColumns=all_orthologs.columns)
    
    # Create Phylogenetic Profiling matrix
    pp = PhylogeneticProfiling(orthologs_swiss_prot_only, onSpecies=["9606"])
    pp_matrix = pp.computePresenceAbscenseMatrix()
    #pp_matrix.to_csv("drive/MyDrive/TFG/orthologs_counts_matrix_v2.tsv", sep="\t", index=True, header=True)
    
    # Assign GO terms
    goa = GeneOntology(obo_file_path="./go.obo", gaf_file_path="./goa_uniprot_qfo.gaf.gz", hasHeader=False)
    goa.filterByAspects([GO_Aspects.BiologicalProcess.value])
    goa.filterByEvidenceCodes(GO_EvidenceCodes.Experimental.value+GO_EvidenceCodes.AuthorStatements.value+["ISS", "RCA", "IC"])
    proteins2GOterms = goa.assignGOterms(pp_matrix.index, include_parents=True)
    pp_matrix["GO_IDs"] = pp_matrix.index.map(lambda x: proteins2GOterms.get(x, pd.NA))
    #pdb.set_trace()
    print(pp_matrix) 
    pp_matrix.to_csv("drive/MyDrive/TFG/phylogenetic_profile_matrix_pres-absc_v2.tab", sep="\t", header=True, index=True)

    # TODO: train_model()

if __name__ == "__main__":
    main()
      
