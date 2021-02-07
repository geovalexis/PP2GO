import os
import logging
from typing import Dict
import pandas as pd
import numpy as np
from enum import Enum, unique
import pdb
import sys


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

    def __init__(self, gaf_file_path: os.path, hasHeader = False):
        self.go_annotations = pd.read_table(gaf_file_path, 
                             header=None, 
                             names=["DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID", "DB:Reference", "Evidence Code", "With (or) From", "Aspect", "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type", "Taxon and Interacting taxon", "Date","Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID"],
                             dtype="string",
                             skiprows=12 if hasHeader else 0, 
                             compression="gzip")
    
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

    def assignGOterms(self, uniprotids: list) -> dict:
        uniprotids = set(uniprotids)
        logging.info(f"Assigning GO terms to {len(uniprotids)} proteins...")
        goa_filteredByUniprotids = self.go_annotations[self.go_annotations["DB_Object_ID"].isin(uniprotids)]
        uniprotids2GOterms = goa_filteredByUniprotids.groupby("DB_Object_ID")["GO_ID"].apply(np.unique).to_dict()
        if len(uniprotids2GOterms) != len(uniprotids):
            logging.warning(f"{len(uniprotids)-len(uniprotids2GOterms)} proteins doesn't have any match.")
        return uniprotids2GOterms

if __name__ == "__main__":
    goa_test = GeneOntology("goa_uniprot_qfo.gaf.gz", hasHeader=False)
    goa_test.filterByAspects([GO_Aspects.BiologicalProcess.value]).filterByEvidenceCodes(GO_EvidenceCodes.Experimental.value+GO_EvidenceCodes.AuthorStatements.value)
    print(goa_test.assignGOterms(["A5PKW4", "O43292"]))
