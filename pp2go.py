#!/usr/bin/env python3

from typing import Dict, List
import pandas as pd
import numpy as np
import os
import logging
import argparse

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

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="PP2GO pipeline", epilog="Enjoy!")
    parser.add_argument("--orthologs", required=True, type=str, help="Orthologs input.")
    parser.add_argument("--gaf-file", required=True, type=str, help="Gene Ontology annotation file.")
    parser.add_argument("--go-obo-file", required=False, default="./data/go.obo", type=str, help="Core ontology in OBO format.")
    parser.add_argument("--output", required=True, type=str, help="Output file.")
    parser.add_argument("-v","--verbose", required=False, default=False, action="store_true", help="Verbose logging.")
    #parser.add_argument("--go-aspects", required=False, type=str.upper, help="GO aspect/ontology", choices=["P", "C", "F"]) #TODO: allow multiple 
    #TODO: add SPECIES, EVIDENCE CODES, Pres-Abs/Counts arguments
    return parser.parse_args()

def main():

    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] -- %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stdout
    )    
    
    # Read and filter orthologs dataset
    all_orthologs = pd.read_table(args.orthologs, names=["uniprotid1", "uniprotid2"], dtype="string")
    orthologs_swiss_prot_only = filterBySwissProt(all_orthologs, onColumns=all_orthologs.columns)
    
    # Compute Phylogenetic Profiling matrix
    pp = PhylogeneticProfiling(orthologs_swiss_prot_only, onSpecies=["9606"])
    pp_matrix = pp.computeCountsMatrix()
    #pp_matrix.to_csv("drive/MyDrive/TFG/orthologs_counts_matrix_v2.tsv", sep="\t", index=True, header=True)
    
    # Assign GO terms
    goa = GeneOntology(obo_file_path=args.go_obo_file, gaf_file_path=args.gaf_file, hasHeader=False)
    goa.filterByAspects([GO_Aspects.BiologicalProcess.value])
    goa.filterByEvidenceCodes(GO_EvidenceCodes.Experimental.value+GO_EvidenceCodes.AuthorStatements.value+["ISS", "RCA", "IC"])
    proteins2GOterms = goa.assignGOterms(pp_matrix.index, include_parents=True)
    pp_matrix["GO_IDs"] = pp_matrix.index.map(lambda x: proteins2GOterms.get(x, np.array([])))
    #pdb.set_trace()

    # Save matrix and GO terms
    pp_matrix["GO_IDs"] = pp_matrix["GO_IDs"].apply(lambda x: ",".join(x))  # Before saving the dataframe we must reformat the lists
    pp_matrix.to_csv(args.output, sep="\t", header=True, index=True)
    print(pp_matrix) 
    
    # TODO: train_model()
    

if __name__ == "__main__":
    main()
      
#python pp2go.py --orthologs drive/MyDrive/TFG/QfO_input.tsv --gaf-file ./data/goa_uniprot_qfo.gaf.gz --output drive/MyDrive/TFG/phylogenetic_profile_matrix_counts_v2.tab -v