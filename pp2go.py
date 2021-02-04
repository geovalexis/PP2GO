#!/usr/bin/env python3

from typing import Dict, List
import pandas as pd
import os

import pdb

__author__ = "Geovanny Risco"
__email__ = "geovanny.risco@bsc.es"
__version__ = "0.1"



__DESCRIPTION__ = ""


def getUniprot2Taxid_NCBI(uniprotIDs: set, prot_accession2taxid: os.path) -> Dict:
    """ Translate a set of genes to its corresponding taxIDs given their uniprotKB accession numbers.
    Args:
        uniprotIDs (set): set of genes to translate
        prot_accession2taxid (os.path): path to the prot.accession2taxid gzip compressed file 
                                        from the official NCBI ftp server.
                                        See ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/ 
    """    
    uniprot2taxid = pd.read_table(prot_accession2taxid,
                            compression="gzip", 
                            names=["accession", "accession.version", "taxid", "gi"],
                            usecols=["accession", "taxid"],
                            dtype={"accession":"string", "taxid": "int32"})
    uniprot2taxid = uniprot2taxid[uniprot2taxid["accession"].isin(uniprotIDs)] 
    # Check if all the uniprotIDs have a corresponding taxID
    if (uniprot2taxid["accession"].nunique() != len(uniprotIDs)):
        print(f"The tax ID for following {len(uniprotIDs)-uniprot2taxid['accession'].nunique()} uniprotKB accession numbers couldn't be found:")
        #print(set.intersection(uniprotIDs, set(uniprot2taxid["accession"].unique())))
    return uniprot2taxid.set_index("accession")["taxid"].to_dict()

def getUniprot2Taxid_Uniprot(uniprotIDs: set, idmapping: os.path) -> Dict:
    """ Translate a set of genes to its corresponding taxIDs given their uniprotKB accession numbers.
    Args:
        uniprotIDs (set): set of genes to translate
        idmapping (os.path): path to the idmapping_selected.tab gzip compressed file 
                             from the official Uniprot ftp server.
                             See ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/idmapping/ 
    """    
    uniprot2taxid = pd.read_table(idmapping,
                            compression="gzip", 
                            names=["UniprotKB-AC", "GO", "NCBI-taxon"],
                            usecols=["UniprotKB-AC", "NCBI-taxon"],
                            dtype={"UniprotKB-AC":"string", "NCBI-taxon": "int32"})
    uniprot2taxid = uniprot2taxid[uniprot2taxid["UniprotKB-AC"].isin(uniprotIDs)] 
    # Check if all the uniprotIDs have a corresponding taxID
    if (uniprot2taxid["UniprotKB-AC"].nunique() != len(uniprotIDs)):
        print(f"The tax ID for following {len(uniprotIDs)-uniprot2taxid['UniprotKB-AC'].nunique()} uniprotKB accession numbers couldn't be found:")
        #print(set.intersection(uniprotIDs, set(uniprot2taxid["UniprotKB-AC"].unique())))
    return uniprot2taxid.set_index("UniprotKB-AC")["NCBI-taxon"].to_dict()


def assignTaxIDs(df: pd.DataFrame, onColumns: List[str], mapping_file: os.path):
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
    uniproid2taxid = getUniprot2Taxid_Uniprot(set(all_uniprotids), mapping_file)
    for column in onColumns:
        df[f"{column}_taxID"] = df[column].apply(lambda x:  uniproid2taxid.get(x, pd.NA)).astype(pd.Int64Dtype())
    return df

def main():
    mtp_orthologs = "drive/MyDrive/TFG/QfO_input.tsv"
    idmapping_qfp_subset = "idmapping_selected_qfo_subset.tab.gz"

    orthologs = pd.read_table(mtp_orthologs, names=["uniprotid1", "uniprotid2"], dtype={"uniprotid1":"string", "uniprotid2":"string"})
    orthologs_withTaxIDs = assignTaxIDs(orthologs, onColumns=["uniprotid1", "uniprotid2"], mapping_file=idmapping_qfp_subset)
    #orthologs_withTaxIDs = orthologs_withTaxIDs[["uniprotid1",  "uniprotid1_taxID", "uniprotid2", "uniprotid2_taxID"]] #Reorder columns
    print(orthologs_withTaxIDs.head())
    
    human_orthologs = ""
    genes2GOterms = None
    # TODO: train_model()
    pass

if __name__ == "__main__":
    main()
      
