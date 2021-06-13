import os
import pdb
import logging
import requests
from enum import Enum
import pandas as pd
import numpy as np

from goatools import obo_parser

__author__ = "Geovanny Risco"
__email__ = "geovanny.risco@bsc.es"
__version__ = "0.1"

logger = logging.getLogger(__name__)


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

    max_level = 16

    def __init__(self, gaf_file_path: os.path, obo_file_path: os.path = None):
        self.go_annotations = pd.read_table(gaf_file_path,                             
                             names=["DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID", "DB:Reference", "Evidence Code", "With (or) From", "Aspect", "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type", "Taxon and Interacting taxon", "Date","Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID"],
                             header=None, 
                             dtype="string",
                             comment="!",
                             compression="gzip")
        self.go = obo_parser.GODag(obo_file_path if obo_file_path else self.retrieveOBOFile(), optional_attrs="relationship") 
    
    @staticmethod
    def retrieveOBOFile():
        go_obo_url = 'http://purl.obolibrary.org/obo/go/go.obo'
        data_folder = os.getcwd() + '/data'
        go_obo_filepath = data_folder+'/go.obo'
        if not os.path.isdir(data_folder):
            os.mkdir(data_folder)
        if not os.path.isfile(go_obo_filepath):
            logger.info("go.obo not found. Downloading last Gene Ontology in OBO format...")
            response = requests.get(go_obo_url)
            with open(go_obo_filepath, "wb") as f:
                f.write(response.content)
        return go_obo_filepath
    
    @property
    def goaDataFrame(self):
        return self.go_annotations

    def filterOutByAspects(self, aspects: list):
        logger.info("Filtering GO terms by aspects...")
        self.go_annotations = self.go_annotations[self.go_annotations["Aspect"].isin(aspects)]
        return self

    def filterOutByEvidenceCodes(self, evidence_codes: list):
        logger.info("Filtering GO terms by evidence codes...")
        self.go_annotations = self.go_annotations[self.go_annotations["Evidence Code"].isin(evidence_codes)]
        return self

    def filterOutByQualifier(self, qualifier: str):
        self.go_annotations = self.go_annotations[self.go_annotations["Qualifier"].fillna("")!=qualifier]
        return self
    
    def setGOtermAsRoot(self, go_term: str):
        logger.info(f"Setting {go_term} as root... All GO terms above it will be filtered out.")
        go_term_children = self.go.get(go_term).get_all_lower() #NOTE:  get_all_children() only get those childrem with is "is_a" relationship (not "part_of", "occurs_in", etc.)
                                                                #       Use .get_all_lower() to get all descendants, included all type of relationships GO terms.
        #go_term_children.add(go_term) #Also add the root GO term as posible candidate
        self.go_annotations = self.go_annotations[self.go_annotations["GO_ID"].isin(go_term_children)]
        return self

    def assignGOterms(self, uniprotids: list, include_parents: bool = False, min_level: int = None, max_level: int = None) -> dict:
        logger.info(f"Assigning GO terms...")
        logger.info(f"{len(uniprotids)} proteins in total...")
        goa_filteredByUniprotids = self.go_annotations[self.go_annotations["DB_Object_ID"].isin(uniprotids)]
        uniprotids2GOterms = goa_filteredByUniprotids.groupby("DB_Object_ID")["GO_ID"].apply(np.unique).to_dict()
        logger.debug(f"{len(uniprotids)-len(uniprotids2GOterms)} proteins were filtered out or were not annotated.")
        # https://nbviewer.jupyter.org/urls/dessimozlab.github.io/go-handbook/GO%20Tutorial%20in%20Python%20-%20Solutions.ipynb
        if include_parents:
            uniprotids2GOterms = self.concatenateParents(uniprotids2GOterms)
        if min_level:
            uniprotids2GOterms = {uniprotid:self.CutOffByMinLevel(go_terms, min_level) for uniprotid, go_terms in uniprotids2GOterms.items()}
        if max_level:
            uniprotids2GOterms = {uniprotid:self.CutOffByMaxLevel(go_terms, max_level) for uniprotid, go_terms in uniprotids2GOterms.items()}

        return uniprotids2GOterms
    
    def findUnannotatedProteins(self):
        raise NotImplementedError()
    
    def CutOffByMaxLevel(self,go_terms: list,  max_level: int) -> list:
        filtered_go_terms = []
        for go_term in go_terms:
            if self.getGOtermLevel(go_term) <= max_level:
                filtered_go_terms.append(go_term)
        return filtered_go_terms
    
    def CutOffByMinLevel(self, go_terms: list, min_level: int) -> list:
        filtered_go_terms = []
        for go_term in go_terms:
            if self.getGOtermLevel(go_term) >= min_level:
                filtered_go_terms.append(go_term)
        return filtered_go_terms

    def getGOtermLevel(self, go_term: str) -> int:
        go_term_record = self.go.get(go_term)
        level = go_term_record.level if go_term_record else self.max_level #NOTE: if no record found, we return the maximum level. REVISE THIS.
        return level
    
    def retrieveParents(self, go_terms: list) -> dict:
        go_terms2parents = {}
        for go_term in set(go_terms):
            go_term_record = self.go.get(go_term)
            go_terms2parents[go_term] = set()
            if go_term_record:
                go_terms2parents[go_term] = go_term_record.get_all_parents() # The built-in function go_term_record.get_all_parents() does not cut-off by level or depth
            else:
                logger.warning(f"Any record has been found for {go_term} in the current OBO file. No parents will be assigned for this GO term.")
        return go_terms2parents

    def concatenateParents(self, uniprotids2GOterms: dict) -> dict:
        all_go_terms = np.hstack(list(uniprotids2GOterms.values())) #Concatenate all GO terms 
        go_terms2parents = self.retrieveParents(all_go_terms)
        unprotids2GOterms_withParents = {}
        for uniprotid, go_terms in uniprotids2GOterms.items():
            go_terms_parents = set()
            for go_term in go_terms:
                go_terms_parents.update(go_terms2parents[go_term])
            unprotids2GOterms_withParents[uniprotid] = np.unique([*go_terms,*go_terms_parents])
        return unprotids2GOterms_withParents


if __name__ == "__main__":
    goa_test = GeneOntology(gaf_file_path="./data/goa_uniprot_qfo.gaf.gz")
    goa_test.filterOutByAspects([GO_Aspects.BiologicalProcess.value]).filterOutByEvidenceCodes(GO_EvidenceCodes.Experimental.value+GO_EvidenceCodes.AuthorStatements.value)
    goa_test.filterOutByQualifier("NOT")
    #goa_test.setGOtermAsRoot("GO:0008150") # GO:0008150 corresponds to the actual root of the Biological Process ontology, which wouldn't affect much if we have already filtered by biological process
    proteins2GOterms = goa_test.assignGOterms(["Q13061", "Q8IZF6", "P08183"], include_parents=True, max_level=1)
    print(proteins2GOterms)
