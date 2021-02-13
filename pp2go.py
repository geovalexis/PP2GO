#!/usr/bin/env python3
from typing import Dict, List
import os
import logging
import argparse
import pdb
import sys

import pandas as pd
import numpy as np

from phylogenetic_profiling import PhylogeneticProfiling
from gene_ontology import GeneOntology, GO_Aspects, GO_EvidenceCodes
from machine_learning import runML

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
    parser.add_argument("--go-obo-file", required=False, type=str, default="./data/go.obo", help="Core ontology in OBO format.")
    parser.add_argument("--pp-matrix", required=False, type=str, default="", help="Name of the Phylogenetic Profiling Matrix if wants to be saved.")
    parser.add_argument("--proteome-species", nargs="*", type=int, default=[9606], help="Space separated list of species whose proteins will be used for the Phylogenetic Profiling Matrix. Human proteome will be taken by default.")
    parser.add_argument("--reference-species", nargs="*", type=int, default=[], help="Space separated list of reference organisms on which the orthologs will be searched for. By default all available will be taken.")
    parser.add_argument("--go-aspects", nargs="*", type=str, default=["P"], choices=["P", "C", "F"], help="GO aspect/ontology. By default only Biological Process will be taken.") 
    parser.add_argument("--ml-results", type=str, required=False, default="", help="Filename for the Machine Learning models assessment results.")
    parser.add_argument("-v","--verbose", required=False, default=False, action="store_true", help="Verbose logging.")
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
    pp = PhylogeneticProfiling(orthologs_swiss_prot_only, onSpecies=args.proteome_species, reference_species=args.reference_species)
    pp_matrix = pp.computeCountsMatrix() #TODO: include Pres-Abs/Counts as arguments
    
    # Assign GO terms
    goa = GeneOntology(obo_file_path=args.go_obo_file, gaf_file_path=args.gaf_file, hasHeader=False) #TODO: support not filtered gaf files
    goa.filterByAspects(args.go_aspects)
    goa.filterByEvidenceCodes(GO_EvidenceCodes.Experimental.value+GO_EvidenceCodes.AuthorStatements.value+["ISS", "RCA", "IC"]) #TODO: add Evidence Code as arguments
    proteins2GOterms = goa.assignGOterms(pp_matrix.index, include_parents=True) #TODO: add include_parents as argument
    pp_matrix["GO_IDs"] = pp_matrix.index.map(lambda x: proteins2GOterms.get(x, np.array([])))
    logging.info(f"Profiling matrix with GO terms...\n{pp_matrix}")
    
    # Performe Machine Learning algorithm (if apply)
    if args.ml_results:
        runML(pp_matrix, min=100, max=1000, results_file=args.ml_results) 

    # Save PP matrix (if apply)
    if args.pp_matrix:
        pp_matrix["GO_IDs"] = pp_matrix["GO_IDs"].apply(lambda x: ",".join(x))  # Before saving the dataframe we must reformat the lists
        pp_matrix.to_csv(args.pp_matrix, sep="\t", header=True, index=True)


if __name__ == "__main__":
    main()


#### TESTS ####
#   
#python pp2go.py --orthologs drive/MyDrive/TFG/QfO_input.tsv --gaf-file ./data/goa_uniprot_qfo.gaf.gz --pp-matrix drive/MyDrive/TFG/results/MTP_last-pp_matrix_counts.tab --ml-results drive/MyDrive/TFG/results/MTP_last-counts-ML_assesment.tab -v
#
#python pp2go.py --orthologs PANTHER_14.1_all-20190603-2359.336.rels.raw --gaf-file ./data/goa_uniprot_qfo.gaf.gz --pp-matrix drive/MyDrive/TFG/results/PANTHER_14.1_all-pp_matrix_counts.tab --ml-results drive/MyDrive/TFG/results/PANTHER_14.1_all-counts-ML_assesment.tab -v