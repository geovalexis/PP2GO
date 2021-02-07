import os
import pdb
import logging
from typing import Dict
from enum import Enum

import pandas as pd
import numpy as np

from goatools import obo_parser





__author__ = "Geovanny Risco"
__email__ = "geovanny.risco@bsc.es"
__version__ = "0.1"


class GO_Aspects(Enum):
    BiologicalProcess = "P"
    CellularComponent = "C"
    MolecularFunction = "F"

class GO_EvidenceCodes(Enum):
    """ See http://geneontology.org/docs/guide-go-evidence-codes/ """
    Experimental = ["EXP", "IMP", "HMP", "IGI", "HGI", "IPI", "IEP", "HEP", "IDA", "HDA"]
    PhylogeneticallyInferred = ["IBA", "IBD", "IKR", "IRD"]
    ComputionalAnalysis = ["ISS", "ISO", "ISA", "ISM", "IGC", "RCA"]
    AuthorStatements = ["TAS", "NAS"]
    CuratorStatements = ["IC", "ND"]
    ElectronicAnnotations = ["IEA"]

class GeneOntology():

    def __init__(self, obo_file_path: os.path, gaf_file_path: os.path, hasHeader: bool):
        self.go_annotations = pd.read_table(gaf_file_path, 
                             header=None, 
                             names=["DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID", "DB:Reference", "Evidence Code", "With (or) From", "Aspect", "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type", "Taxon and Interacting taxon", "Date","Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID"],
                             dtype="string",
                             skiprows=12 if hasHeader else 0, 
                             compression="gzip")
        self.go = obo_parser.GODag(obo_file_path)
    
    @property
    def goaDataFrame(self):
        return self.go_annotations

    def filterByAspects(self, aspects: list):
        logging.info("Filtering GO terms by aspects...")
        self.go_annotations = self.go_annotations[self.go_annotations["Aspect"].isin(aspects)]
        return self

    def filterByEvidenceCodes(self, evidence_codes: list):
        logging.info("Filtering GO terms by evidence codes...")
        self.go_annotations = self.go_annotations[self.go_annotations["Evidence Code"].isin(evidence_codes)]
        return self

    def assignGOterms(self, uniprotids: list, include_parents: bool = False) -> dict:
        uniprotids = set(uniprotids)
        logging.info(f"Assigning GO terms to {len(uniprotids)} proteins...")
        goa_filteredByUniprotids = self.go_annotations[self.go_annotations["DB_Object_ID"].isin(uniprotids)]
        uniprotids2GOterms = goa_filteredByUniprotids.groupby("DB_Object_ID")["GO_ID"].apply(np.unique).to_dict()
        if len(uniprotids2GOterms) != len(uniprotids):
            logging.warning(f"{len(uniprotids)-len(uniprotids2GOterms)} proteins doesn't have any match.")
        if include_parents:
            uniprotids2GOterms = self.concatenateParents(uniprotids2GOterms)
        return uniprotids2GOterms
    
    def retrieveParents(self, go_terms: list) -> dict:
        go_terms2parents = {}
        for go_term in set(go_terms):
            go_term_record = self.go.get(go_term)
            if go_term_record:
                go_terms2parents[go_term] = go_term_record.get_all_parents()
            else:
                logging.warning(f"Any record has been found for {go_term} in the current OBO file. No parents will be assigned for this GO term.")
                go_terms2parents[go_term] = set()
        return go_terms2parents

    def concatenateParents(self, uniprotids2GOterms: dict) -> dict:
        all_go_terms = np.hstack(list(uniprotids2GOterms.values())) #Concatenate all GO terms 
        go_terms2parents = self.retrieveParents(all_go_terms)
        unprotids2GOterms_withParents = {}
        for uniprotid, go_terms in uniprotids2GOterms.items():
            go_terms_parents = set()
            for go_term in go_terms:
                go_terms_parents.update(go_terms2parents[go_term])
            unprotids2GOterms_withParents[uniprotid] = np.hstack([go_terms,*go_terms_parents])
        return unprotids2GOterms_withParents

if __name__ == "__main__":
    goa_test = GeneOntology(obo_file_path="./go.obo", gaf_file_path="./goa_uniprot_qfo.gaf.gz", hasHeader=False)
    goa_test.filterByAspects([GO_Aspects.BiologicalProcess.value]).filterByEvidenceCodes(GO_EvidenceCodes.Experimental.value+GO_EvidenceCodes.AuthorStatements.value)
    proteins2GOterms = goa_test.assignGOterms(["A5PKW4", "O43292", "P08183"], include_parents=True)
    print(proteins2GOterms)
