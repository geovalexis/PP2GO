import os
import logging
from typing import Dict
import pandas as pd

import aiohttp
import asyncio
import pdb
import sys


__author__ = "Geovanny Risco"
__email__ = "geovanny.risco@bsc.es"
__version__ = "0.1"



__DESCRIPTION__ = ""

class PhylogeneticProfiling():

    def __init__(self, orthologs: pd.DataFrame, onSpecies: list, reference_species: list = None):
        self._orthologs = self.mapTaxIDs(orthologs, onColumns=orthologs.columns, dropUnmatched=True)
        self._onSpecies = onSpecies
        self._reference_species = reference_species

    # NOTE: we cannot define multiple constructors in Python
    # def __init__(self, orthologs: pd.DataFrame, onProteins: list, reference_species: list, mapping_file: os.path):
    #     #TODO: we might need to compute a phylogentic profile given a list of proteins instead of a list of species
    #     raise NotImplementedError()

    
    @property
    def orthologsDataFrame(self):
        return self._orthologs
    
    @property
    def referenceSpecies(self):
        return self._reference_species
    
    @property
    def onSpecies(self):
        return self._onSpecies

    @staticmethod
    def mapTaxIDs(df: pd.DataFrame, onColumns: list, dropUnmatched: bool = True):
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
        #uniproid2taxid = loop.run_until_complete(TaxaMapping.mapUniprotIDs2taxIDs_EBIRequest_multiprocessing(list(set(all_uniprotids))))
        #uniproid2taxid = TaxaMapping.mapUniprot2Taxid_NCBI(set(all_uniprotids))
        uniproid2taxid = TaxaMapping.mapUniprot2Taxid_Uniprot(set(all_uniprotids))
        for column in onColumns:
            df[f"{column}_taxID"] = df[column].apply(lambda x:  uniproid2taxid.get(x, pd.NA)).astype(pd.Int64Dtype())
        if dropUnmatched: 
            df.dropna()
        return df



class TaxaMapping():

    @staticmethod
    def mapUniprot2Taxid_NCBI(uniprotIDs: set, prot_accession2taxid: os.path = "prot.accession2QfOtaxid.gz") -> Dict:
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

    
    @staticmethod
    def mapUniprot2Taxid_Uniprot(uniprotIDs: set, idmapping: os.path = "idmapping_selected_qfo_subset.tab.gz") -> Dict: 
        """ Translate a set of genes to its corresponding taxIDs given their uniprotKB accession numbers.
        Args:
            uniprotIDs (set): set of genes to translate
            idmapping (os.path): path to the idmapping_selected.tab gzip compressed file 
                                from the official Uniprot ftp server.
                                See ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/idmapping/ 
        """    
        uniprot2taxid = pd.read_table(idmapping,
                                compression="gzip", 
                                names=["UniprotKB-AC", "GO", "NCBI-taxon"],  #TODO: use the original file with all the columns
                                usecols=["UniprotKB-AC", "NCBI-taxon"],
                                dtype={"UniprotKB-AC":"string", "NCBI-taxon": "int32"})
        uniprot2taxid = uniprot2taxid[uniprot2taxid["UniprotKB-AC"].isin(uniprotIDs)] 
        # Check if all the uniprotIDs have a corresponding taxID
        if (uniprot2taxid["UniprotKB-AC"].nunique() != len(uniprotIDs)):
            logging.warning(f"The tax IDs for {len(uniprotIDs)-uniprot2taxid['UniprotKB-AC'].nunique()} uniprotKB accession numbers couldn't be found.")
            logging.debug(set.difference(uniprotIDs, set(uniprot2taxid["UniprotKB-AC"].unique())))
        return uniprot2taxid.set_index("UniprotKB-AC")["NCBI-taxon"].to_dict()

    @staticmethod
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
                
    @classmethod
    async def mapUniprotIDs2taxIDs_EBIRequest_multiprocessing(cls, uniprotIDs, chunk=200) -> Dict: # EBI limits requests to 200 requests/second/user
        result = {}
        async with aiohttp.ClientSession() as session:
            if len(uniprotIDs)>chunk:
                for i in range(0, len(uniprotIDs), chunk):
                    step = i+chunk if (i+chunk) < len(uniprotIDs) else len(uniprotIDs)
                    logging.info(f"{i} proteins already processed. Processing next batch...")
                    res_batch = await asyncio.gather(*[cls.map_uniprotIDs2taxIDs_EBIRequest(uniprotID, session) for uniprotID in uniprotIDs[i:step]])
                    for j in res_batch:
                        result.update(j)
                    #sleep(1)
            else:
                res_batch = await asyncio.gather(*[cls.map_uniprotIDs2taxIDs_EBIRequest(uniprotID, session) for uniprotID in uniprotIDs])
                for j in res_batch:
                    result.update(j)
        if (len(uniprotIDs) != len(result)):
            logging.warning(f"The tax IDs for {len(uniprotIDs)-len(result.keys())} uniprotKB accession numbers couldn't be found.")
            logging.debug(set.difference(set(uniprotIDs), result.keys()))
        return result    
