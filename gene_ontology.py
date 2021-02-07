import os
import logging
from typing import Dict
import pandas as pd
import numpy as np
from enum import Enum
import pdb
import sys


__author__ = "Geovanny Risco"
__email__ = "geovanny.risco@bsc.es"
__version__ = "0.1"


class GOAspects(Enum):
    BiologicalProcess = "P"
    CellularComponent = "C"
    MolecularFunction = "F"

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

    def filterByAspect(self, aspect: GOAspects):
        self.go_annotations = self.go_annotations[self.go_annotations["Aspect"]==aspect.value]
        return self




if __name__ == "__main__":
    test = GeneOntology("goa_uniprot_qfo.gaf.gz", hasHeader=False)
    test.filterByAspect(GOAspects.BiologicalProcess)
    
    print(test.goaDataFrame)
      