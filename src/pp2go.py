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
from helpers.helper_functions import downloadSwissProtIds

__author__ = "Geovanny Risco"
__email__ = "geovanny.risco@bsc.es"
__version__ = "0.1"
__DESCRIPTION__ = ""


IDMAPPING_FILEPATH = "./data/idmapping_selected_qfo_subset.tab.gz"
GAF_FILEPATH = "./data/goa_uniprot_qfo.gaf.gz"


def intersectSwissProt(uniprotIDS: set, swiss_prot_ids: os.path):
    with open(swiss_prot_ids, "r") as input:
        sp_ids = input.read().split()
    sp_ids = list(map(lambda line: line.split(","), sp_ids))
    sp_ids_exploded = [item for packed_elements in sp_ids for item in packed_elements]
    return uniprotIDS.intersection(set(sp_ids_exploded))

def filterBySwissProt(df: pd.DataFrame, onColumns: List[str]):
    data_folder = os.getcwd() + '/data'
    sp_ids_filepath = data_folder+'/swiss_prot_ids.txt'
    if not os.path.isdir(data_folder):
        os.mkdir(data_folder)
    if not os.path.isfile(sp_ids_filepath):
        logging.info("SwissProt Identifiers list not found. Downloading the last SwissProt dataset from Uniprot...")
        sp_ids = downloadSwissProtIds()
        with open(sp_ids_filepath, "w") as out_sp:
            out_sp.write("\n".join(sp_ids))
    else:
        with open(sp_ids_filepath, "r") as swiss_prot_ids_file:
            sp_ids = swiss_prot_ids_file.read().splitlines()

    for column in onColumns:
        df = df[df[column].isin(sp_ids)]
    return df   


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="PP2GO pipeline", epilog="Enjoy!")
    parser.add_argument("--gaf-file", required=True if not os.path.isfile(GAF_FILEPATH) else False, type=str, default=GAF_FILEPATH, help="Gene Ontology annotation file.")
    parser.add_argument("--idmapping_file", required=True if not os.path.isfile(IDMAPPING_FILEPATH) else False, type=str, default=IDMAPPING_FILEPATH, help="Identifiers mapping file from Uniprot.")
    parser.add_argument("--orthologs", required=True, type=str, help="Orthologs input.")
    parser.add_argument("--filter-by-sp", required=False, default=False, action="store_true", help="Select if wants to filter by only Swiss Prot proteins")
    parser.add_argument("--pres-abs", required=False, default=False, action="store_true", help="Compute presence/abscense PP Matrix instead of counts PP matrix.")
    parser.add_argument("--proteome-species", nargs="*", type=int, default=[9606], help="Space separated list of species whose proteins will be used for the Phylogenetic Profiling Matrix. Human proteome will be taken by default.")
    parser.add_argument("--proteome", type=argparse.FileType("rt"), required=False, help="File that contains a line-separated list of proteins to be take as input proteome.")
    parser.add_argument("--reference-species", nargs="*", type=int, default=[], help="Space separated list of reference organisms on which the orthologs will be searched for. By default all available will be taken.")
    parser.add_argument("--pp-matrix", required=False, type=str, default="", help="Name of the Phylogenetic Profiling Matrix if wants to be saved.")
    parser.add_argument("--go-aspects", nargs="*", type=str, default=["P"], choices=["P", "C", "F"], help="GO aspect/ontology. By default only Biological Process will be taken.") 
    parser.add_argument("--set-as-root", type=str, required=False, default=None, help="Set a given GO term as root, so only their children can be assigned")
    parser.add_argument("--include-go-parents", required=False, default=False, action="store_true", help="Include all the lineage for each GO term assign to a protein.")
    parser.add_argument("--min-level", type=int, required=False, default=None, help="Minimum level until which the parents will be included (if apply)")
    parser.add_argument("--min-depth", type=int, required=False, default=None, help="Minimum depth until which the parents will be included (if apply)")
    parser.add_argument("--min-gos", type=int, required=False, default=None, help="Min number of GO terms' ocurrences,")
    parser.add_argument("--max-gos", type=int, required=False, default=None, help="Max number of GO terms' ocurrences,")
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
    orthologs = pd.read_table(args.orthologs, names=["uniprotid1", "uniprotid2"], dtype="string")
    if args.filter_by_sp: orthologs = filterBySwissProt(orthologs, onColumns=orthologs.columns) 

    # Compute Phylogenetic Profiling matrix
    pp = PhylogeneticProfiling(idmapping_file=args.idmapping_file, 
                                orthologs=orthologs, 
                                reference_species=args.reference_species, 
                                onProteins=args.proteome.read().splitlines() if args.proteome else [], 
                                onSpecies=args.proteome_species)
    pp_matrix = pp.computeCountsMatrix() if not args.pres_abs else pp.computePresenceAbscenseMatrix() 
    
    # Assign GO terms
    goa = GeneOntology(gaf_file_path=args.gaf_file)
    goa.filterByAspects(args.go_aspects)
    goa.filterByQualifier("NOT") #TODO: add as argument
    goa.filterByEvidenceCodes(GO_EvidenceCodes.Experimental.value+GO_EvidenceCodes.AuthorStatements.value+["ISS", "RCA", "IC"]) #TODO: add Evidence Code as arguments
    if args.set_as_root: goa.setGOtermAsRoot(args.set_as_root)
    # TODO: save proteins that are not really annotated -> in this case an empty list does not mean that proteins is not annotated because it is not the original GAF dataset (it has been filtered)
    proteins2GOterms = goa.assignGOterms(pp_matrix.index.get_level_values(1), include_parents=args.include_go_parents, min_level=args.min_level, min_depth=args.min_depth) 
    pp_matrix["GO_IDs"] = pp_matrix.index.map(lambda x: proteins2GOterms.get(x[1], np.array([])))
    logging.info(f"Profiling matrix with GO terms...\n{pp_matrix}")
    
    
    try:
        # Performe Machine Learning algorithm (if apply)
        runML(pp_matrix, min=args.min_gos, max=args.max_gos, results_file=args.ml_results) 
    except:
        raise
    finally:
        # Save PP matrix (if apply)
        if args.pp_matrix:
            pp_matrix["GO_IDs"] = pp_matrix["GO_IDs"].apply(lambda x: ",".join(x))  # Before saving the dataframe we must reformat the lists
            pp_matrix.to_csv(args.pp_matrix, sep="\t", header=True, index=True)


if __name__ == "__main__":
    main()


#### TESTS ####
#   
#python pp2go.py --orthologs drive/MyDrive/TFG/QfO_input.tsv --include-go-parents --pp-matrix drive/MyDrive/TFG/results/MTP_last-pp_matrix_counts.tab --ml-results drive/MyDrive/TFG/results/MTP_last-counts-ML_assesment.tab --min-gos 100 --min-level 3 -v
#
#python pp2go.py --orthologs PANTHER_14.1_all-20190603-2359.336.rels.raw --include-go-parents--pp-matrix drive/MyDrive/TFG/results/PANTHER_14.1_all-pp_matrix_counts.tab --ml-results drive/MyDrive/TFG/results/PANTHER_14.1_all-counts-ML_assesment.tab --min-gos 100 --min-level 3 -v