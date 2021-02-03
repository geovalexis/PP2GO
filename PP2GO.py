#!/usr/bin/env python3

from typing import Dict
import pandas as pd
import gzip

__author__ = "Geovanny Risco"
__email__ = "geovanny.risco@bsc.es"
__version__ = "0.1"



__DESCRIPTION__ = ""


def getTaxIDs(uniprotIDs: set) -> Dict:
    """ Translate a set of genes to its corresponding taxIDs given their uniprotKB accession numbers.
    Args:
        uniprotIDs (set): set of genes to translate
    """
    uniprot2taxid = pd.read_table("prot.accession2tQfOaxid.gz",
                            compression="gzip", 
                            names=["accession", "accession.version", "taxid", "gi"],
                            usecols=["accession", "taxid"],
                            dtype={"accession":"string", "taxid": "int32"})
    uniprot2taxid = uniprot2taxid[uniprot2taxid["accession"].isin(uniprotIDs)] 
    return uniprot2taxid.set_index("accession")["taxid"].to_dict()

def main():
    orthologs = None
    genes2GOterms = None
    # TODO: train_model()
    pass

if __name__ == "__main__":
    main()
      
