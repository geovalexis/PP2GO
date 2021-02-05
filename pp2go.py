#!/usr/bin/env python3

from typing import Dict, List
import pandas as pd
import os
import logging
import requests
import aiohttp
import asyncio
import pdb
import sys

__author__ = "Geovanny Risco"
__email__ = "geovanny.risco@bsc.es"
__version__ = "0.1"



__DESCRIPTION__ = ""


def mapUniprot2Taxid_NCBI(uniprotIDs: set, prot_accession2taxid: os.path) -> Dict:
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
        logging.warning(f"The tax ID for {len(uniprotIDs)-uniprot2taxid['accession'].nunique()} uniprotKB accession numbers couldn't be found.")
        logging.debug(set.difference(uniprotIDs, set(uniprot2taxid["accession"].unique())))
    return uniprot2taxid.set_index("accession")["taxid"].to_dict()

def mapUniprot2Taxid_Uniprot(uniprotIDs: set, idmapping: os.path) -> Dict:
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
        logging.warning(f"The tax IDs for {len(uniprotIDs)-uniprot2taxid['UniprotKB-AC'].nunique()} uniprotKB accession numbers couldn't be found.")
        logging.debug(set.difference(uniprotIDs, set(uniprot2taxid["UniprotKB-AC"].unique())))
    return uniprot2taxid.set_index("UniprotKB-AC")["NCBI-taxon"].to_dict()

async def map_uniprotIDs2taxIDs_EBIRequest(uniprotID, session):
    # Documentation in https://www.ebi.ac.uk/proteins/api/doc/
    while True:
        try:
            requestURL = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprotID}"
            response = await session.get(requestURL, headers={ "Accept" : "application/json"})
            if response.ok: # status_code == 200
                return {uniprotID:(await response.json())["organism"]["taxonomy"]}
            else:
                logging.debug(f"Protein {uniprotID} raised a {response.status} status code.") # It means that the protein couldn't be found or has been deleted.
                return {}
        except aiohttp.ClientConnectionError as e:
            logging.error(f"Raised a ClientConnectionError: {e.message}")
            asyncio.sleep(0.1)
            logging.info(f"Retrying protein {uniprotID}")
            

async def mapUniprotIDs2taxIDs_EBIRequest_multiprocessing(uniprotIDs, chunk=200) -> Dict: # EBI limits requests to 200 requests/second/user
    result = {}
    async with aiohttp.ClientSession() as session:
        if len(uniprotIDs)>chunk:
            for i in range(0, len(uniprotIDs), chunk):
                step = i+chunk if (i+chunk) < len(uniprotIDs) else len(uniprotIDs)
                logging.info(f"{i} proteins already processed. Processing next batch...")
                res_batch = await asyncio.gather(*[map_uniprotIDs2taxIDs_EBIRequest(uniprotID, session) for uniprotID in uniprotIDs[i:step]])
                for j in res_batch:
                    result.update(j)
                #sleep(1)
        else:
            res_batch = await asyncio.gather(*[map_uniprotIDs2taxIDs_EBIRequest(uniprotID, session) for uniprotID in uniprotIDs])
            for j in res_batch:
                result.update(j)
    if (len(uniprotIDs) != len(result)):
        logging.warning(f"The tax IDs for {len(uniprotIDs)-len(result.keys())} uniprotKB accession numbers couldn't be found.")
        logging.debug(set.difference(set(uniprotIDs), result.keys()))
    return result    

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
    logging.info("Starting taxID assignment...")
    #loop = asyncio.get_event_loop()
    #uniproid2taxid = loop.run_until_complete(mapUniprotIDs2taxIDs_EBIRequest_multiprocessing(list(set(all_uniprotids))))
    uniproid2taxid = mapUniprot2Taxid_Uniprot(set(all_uniprotids), mapping_file)
    for column in onColumns:
        df[f"{column}_taxID"] = df[column].apply(lambda x:  uniproid2taxid.get(x, pd.NA)).astype(pd.Int64Dtype())
    return df


def intersectSwissProt(uniprotIDS: set, swiss_prot_ids: os.path):
    with open(swiss_prot_ids, "r") as input:
        sp_ids = input.read().split()
    sp_ids = list(map(lambda line: line.split(","), sp_ids))
    sp_ids_exploded = [item for packed_elements in sp_ids for item in packed_elements]
    return uniprotIDS.intersection(set(sp_ids_exploded))

def filterBySwissProt(df: pd.DataFrame, onColumns: List[str], swiss_prot_ids: str):
    for column in onColumns:
        swiss_prot_proteins = intersectSwissProt(set(df[column].unique()), swiss_prot_ids)
        df = df[df[column].isin(swiss_prot_proteins)]
    return df

def main():
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] -- %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stdout
    )    

    mtp_orthologs = "drive/MyDrive/TFG/QfO_input.tsv"
    idmapping_qfp_subset = "idmapping_selected_qfo_subset.tab.gz"
    swiss_prot_proteins = "drive/MyDrive/TFG/swiss_prot_ids.csv"

    orthologs = pd.read_table(mtp_orthologs, names=["uniprotid1", "uniprotid2"], dtype={"uniprotid1":"string", "uniprotid2":"string"})
    orthologs_swiss_prot_only = filterBySwissProt(orthologs, onColumns=["uniprotid1", "uniprotid2"], swiss_prot_ids=swiss_prot_proteins)
    orthologs_sp_withTaxIDs = assignTaxIDs(orthologs_swiss_prot_only, onColumns=["uniprotid1", "uniprotid2"], mapping_file=idmapping_qfp_subset)
    orthologs_sp_withTaxIDs = orthologs_sp_withTaxIDs[["uniprotid1",  "uniprotid1_taxID", "uniprotid2", "uniprotid2_taxID"]] #Reorder columns
    #pdb.set_trace()
    print(orthologs_sp_withTaxIDs.head())
    
    human_orthologs = ""
    genes2GOterms = None
    # TODO: train_model()
    pass

if __name__ == "__main__":
    main()
      
