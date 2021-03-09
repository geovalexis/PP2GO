import os
import pdb
import logging
from numpy.lib.function_base import kaiser
import requests
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

    def __init__(self, gaf_file_path: os.path, obo_file_path: os.path = None):
        self.go_annotations = pd.read_table(gaf_file_path,                             
                             names=["DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID", "DB:Reference", "Evidence Code", "With (or) From", "Aspect", "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type", "Taxon and Interacting taxon", "Date","Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID"],
                             header=None, 
                             dtype="string",
                             comment="!",
                             compression="gzip")
        self.go = obo_parser.GODag(obo_file_path if obo_file_path else self.retrieveOBOFile()) 
    
    @staticmethod
    def retrieveOBOFile():
        go_obo_url = 'http://purl.obolibrary.org/obo/go/go.obo'
        data_folder = os.getcwd() + '/data'
        go_obo_filepath = data_folder+'/go.obo'
        if not os.path.isdir(data_folder):
            os.mkdir(data_folder)
        if not os.path.isfile(go_obo_filepath):
            logging.info("go.obo not found. Downloading last Gene Ontology in OBO format...")
            response = requests.get(go_obo_url)
            with open(go_obo_filepath, "wb") as f:
                f.write(response.content)
        return go_obo_filepath
    
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

    def filterByQualifier(self, qualifier: str):
        self.go_annotations = self.go_annotations[self.go_annotations["Qualifier"].fillna("")!=qualifier]
        return self
    
    def setGOtermAsRoot(self, go_term: str):
        logging.info(f"Setting {go_term} as root... All GO terms above it will be filtered out.")
        go_term_children = self.go.get(go_term).get_all_children()
        self.go_annotations = self.go_annotations[self.go_annotations["GO_ID"].isin(go_term_children)]
        return self

    def assignGOterms(self, uniprotids: list, include_parents: bool = False, **kwargs) -> dict:
        logging.info(f"Assigning GO terms...")
        logging.info(f"{len(uniprotids)} proteins in total...")
        goa_filteredByUniprotids = self.go_annotations[self.go_annotations["DB_Object_ID"].isin(uniprotids)]
        uniprotids2GOterms = goa_filteredByUniprotids.groupby("DB_Object_ID")["GO_ID"].apply(np.unique).to_dict()
        #TODO: instead of taking all GO terms assigned, look for the deepest common ancestor between them: 
        # https://nbviewer.jupyter.org/urls/dessimozlab.github.io/go-handbook/GO%20Tutorial%20in%20Python%20-%20Solutions.ipynb
        if len(uniprotids2GOterms) != len(uniprotids):
            logging.warning(f"{len(uniprotids)-len(uniprotids2GOterms)} proteins doesn't have any GO term.")
        if include_parents:
            uniprotids2GOterms = self.concatenateParents(uniprotids2GOterms, **kwargs)
        return uniprotids2GOterms
    
    def retrieveParents(self, go_terms: list, **kwargs) -> dict:
        if kwargs.get("min_level") and kwargs.get("min_depth"):
            raise Exception("You must select only one of the options: either max_level or max_depth, not both")
        go_terms2parents = {}
        for go_term in set(go_terms):
            go_term_record = self.go.get(go_term)
            go_terms2parents[go_term] = set()
            if go_term_record:
                self.findParents(go_term_record, go_terms2parents[go_term], **kwargs) # The built-in function go_term_record.get_all_parents() does not cut-off by level or depth
            else:
                logging.warning(f"Any record has been found for {go_term} in the current OBO file. No parents will be assigned for this GO term.")
        return go_terms2parents
    
    def findParents(self, term_record, parents: set, min_level: int = None, min_depth: int = None):
        if min_level and term_record.level < min_level:
            return
        elif min_depth and term_record.depth < min_depth:
            return
        for parent_term in term_record.parents:
            parents.update({parent_term.item_id})
            self.findParents(parent_term, parents, min_level, min_depth)

    def concatenateParents(self, uniprotids2GOterms: dict, **kwargs) -> dict:
        all_go_terms = np.hstack(list(uniprotids2GOterms.values())) #Concatenate all GO terms 
        go_terms2parents = self.retrieveParents(all_go_terms, **kwargs)
        unprotids2GOterms_withParents = {}
        for uniprotid, go_terms in uniprotids2GOterms.items():
            go_terms_parents = set()
            for go_term in go_terms:
                go_terms_parents.update(go_terms2parents[go_term])
            unprotids2GOterms_withParents[uniprotid] = np.unique([*go_terms,*go_terms_parents])
        return unprotids2GOterms_withParents


if __name__ == "__main__":
    goa_test = GeneOntology(gaf_file_path="./data/goa_uniprot_qfo.gaf.gz")
    goa_test.filterByAspects([GO_Aspects.BiologicalProcess.value]).filterByEvidenceCodes(GO_EvidenceCodes.Experimental.value+GO_EvidenceCodes.AuthorStatements.value)
    goa_test.filterByQualifier("NOT")
    goa_test.setGOtermAsRoot("GO:0008150") # GO:0008150 corresponds to the actual root of the Biological Process ontology, which wouldn't affect much if we have already filtered by biological process
    proteins2GOterms = goa_test.assignGOterms(["Q13061", "Q8IZF6", "P08183"], include_parents=True, min_level=3, min_depth=None)
    print(proteins2GOterms)
