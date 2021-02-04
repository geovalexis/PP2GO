#!/usr/bin/env python3

from typing import Dict, List
import pandas as pd
import os

import pdb

__author__ = "Geovanny Risco"
__email__ = "geovanny.risco@bsc.es"
__version__ = "0.1"



__DESCRIPTION__ = ""

def getUniprot2Taxid(uniprotIDs: set, prot_accession2taxid: os.path) -> Dict:
    """ Translate a set of genes to its corresponding taxIDs given their uniprotKB accession numbers.
    Args:
        uniprotIDs (set): set of genes to translate
        prot_accession2taxid (os.path): path to the prot.accession2taxid gzip compressed file 
                                        from the official NCBI ftp server.
    """    
    uniprot2taxid = pd.read_table(prot_accession2taxid,
                            compression="gzip", 
                            names=["accession", "accession.version", "taxid", "gi"],
                            usecols=["accession", "taxid"],
                            dtype={"accession":"string", "taxid": "int32"})
    uniprot2taxid = uniprot2taxid[uniprot2taxid["accession"].isin(uniprotIDs)] 
    # Check if all the uniprotIDs have a corresponding taxID
    if (uniprot2taxid["accession"].nunique() != len(uniprotIDs)):
        #pdb.set_trace()
        print(f"The tax ID for following {len(uniprotIDs)-uniprot2taxid['accession'].nunique()} uniprotKB accession numbers couldn't be found:")
        #print(set.intersection(uniprotIDs, set(uniprot2taxid["accession"].unique())))
    return uniprot2taxid.set_index("accession")["taxid"].to_dict()


def assignTaxIDs(df: pd.DataFrame, onColumns: List[str], prot_accession2taxid: os.path):
    """[summary]

    Args:
        df (pd.DataFrame): [description]
        onColumns (List[str]): [description]

    Returns:
        df: same dataframe as passed in but with the taxIDs columns added with 
            the following format <onColumn>_taxID
    """
    all_uniprotids = []
    for column in onColumns:
        all_uniprotids.extend(df[column].unique())
    #pdb.set_trace()
    uniproid2taxid = getUniprot2Taxid(set(all_uniprotids), prot_accession2taxid)
    for column in onColumns:
        df[f"{column}_taxID"] = df[column].apply(lambda x:  uniproid2taxid.get(x, pd.NA)).astype(pd.Int64Dtype())
    return df

def main():
    mtp_orthologs = "drive/MyDrive/TFG/QfO_input.tsv"
    qfo_protAccession2taxid = "prot.accession2QfOtaxid.gz"

    orthologs = pd.read_table(mtp_orthologs, names=["uniprotid1", "uniprotid2"], dtype={"uniprotid1":"string", "uniprotid2":"string"})
    orthologs_withTaxIDs = assignTaxIDs(orthologs, onColumns=["uniprotid2"], prot_accession2taxid=qfo_protAccession2taxid)
    #orthologs_withTaxIDs = orthologs_withTaxIDs[["uniprotid1",  "uniprotid1_taxID", "uniprotid2", "uniprotid2_taxID"]] #Reorder columns
    print(orthologs_withTaxIDs.head())
    
    human_orthologs = ""
    genes2GOterms = None
    # TODO: train_model()
    pass

if __name__ == "__main__":
    main()
      
