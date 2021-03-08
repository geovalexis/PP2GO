import os
import logging
from typing import Dict
import pandas as pd
import numpy as np
import aiohttp
import asyncio
import pdb
import sys
from joblib import Parallel, delayed, parallel

__author__ = "Geovanny Risco"
__email__ = "geovanny.risco@bsc.es"
__version__ = "0.1"



__DESCRIPTION__ = ""

class PhylogeneticProfiling():
    """ 
    It takes an orthologs dataset as input in the following format: uniprotid1 | uniprotid2 | uniprotid1_taxID | uniprot2_taxID. 
    Both set of proteins are supposed to be interexchangable, so in this case the first set of uniprot IDs will act as proteome 
    dataset and the second as orthologs dataset. 
    """

    def __init__(self, idmapping_file: os.path, orthologs: pd.DataFrame, reference_species: list = [], onProteins: list = [], onSpecies: list = []):
        self._proteome_column = orthologs.columns[0]
        self._reference_species_column = orthologs.columns[1]
        self._species_availables = None
        self._orthologs = self.mapAndfilterOrthologs(idmapping_file, orthologs, reference_species, onProteins, onSpecies)

    def mapAndfilterOrthologs(self, idmapping_file, orthologs, reference_species, onProteins, onSpecies):
        """
        Map taxIDs to the orthologs and filter dataset by the given reference species and proteins/species
        """
        taxd_orthologs = self.mapTaxIDs(idmapping_file, orthologs, onColumns=[self._proteome_column, self._reference_species_column], dropUnmatched=True) #TODO: Support for OrthoXML format
        if reference_species:
            filtered_orthologs_dataset = self.filterByReferenceSpecies(taxd_orthologs, set(reference_species))      
        else:
            filtered_orthologs_dataset = taxd_orthologs  

        if onProteins:
            filtered_orthologs_dataset = self.filterByProteins(taxd_orthologs, set(onProteins))
        elif onSpecies:
            filtered_orthologs_dataset = self.filterBySpecies(taxd_orthologs,  set(onSpecies))
        #TODO: check that there is no repeated orthologs (bot uniprotId1 and uniprotId2 are the same)
        return filtered_orthologs_dataset

    def filterByReferenceSpecies(self, taxd_orthologs_dataset, reference_species):
        logging.debug("Filtering orthologs dataset by reference species...")
        self._species_availables = list(taxd_orthologs_dataset[f"{self._reference_species_column}_taxID"].unique()) #Save all species available in the dataset
        return taxd_orthologs_dataset[taxd_orthologs_dataset[f"{self._reference_species_column}_taxID"].isin(reference_species)]

    def filterByProteins(self, taxd_orthologs_dataset, onProteins):
        logging.debug("Filtering orthologs dataset by proteins...")
        return taxd_orthologs_dataset[taxd_orthologs_dataset[self._proteome_column].isin(onProteins)]

    def filterBySpecies(self, taxd_orthologs_dataset, onSpecies):
        logging.debug("Filtering orthologs dataset by species...")
        return taxd_orthologs_dataset[taxd_orthologs_dataset[f"{self._proteome_column}_taxID"].isin(onSpecies)]

    @property
    def referenceSpecies(self):
        return list(self._orthologs[f"{self._reference_species_column}_taxID"].unique())
    
    @property
    def onSpecies(self):
        return list(self._orthologs[f"{self._proteome_column}_taxID"].unique())

    @property
    def onProteins(self):
        return list(self._orthologs[self._proteome_column].unique())
    

    @staticmethod
    def mapTaxIDs(idmapping: os.path, df: pd.DataFrame, onColumns: list, dropUnmatched: bool = True):
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
        uniproid2taxid = TaxaMapping.mapUniprot2Taxid_Uniprot(set(all_uniprotids), idmapping)
        for column in onColumns:
            df[f"{column}_taxID"] = df[column].apply(lambda x:  uniproid2taxid.get(x, pd.NA)).astype(pd.Int64Dtype())
        if dropUnmatched: 
            df.dropna(inplace=True)
        return df

    def computePresenceAbscenseMatrix(self):
        return self.computeCountsMatrix().applymap(lambda x: 1 if x >= 1 else 0)

    def computeCountsMatrix(self):
        proteome = list(self._orthologs.drop_duplicates(subset=[self._proteome_column])[self._proteome_column])
        proteome_taxa = list(self._orthologs.drop_duplicates(subset=[self._proteome_column])[f"{self._proteome_column}_taxID"])
        logging.debug(f"Total number of unique proteins: {len(proteome)}")
        logging.info("Computing Phylogenetic Profiling matrix...")
        with Parallel(n_jobs=-1) as run_in_parallel:
            matrix = dict(filter(None, run_in_parallel(delayed(self.searchOrtholog)(tax, self._orthologs, proteome) for tax in self.referenceSpecies)))
        logging.info(f"...found orthologs in {len(matrix)} out of {len(self._species_availables)} taxons.")
        matrix_df = pd.DataFrame(matrix, index=[proteome_taxa,proteome])
        logging.info(f"Final shape of the matrix: {matrix_df.shape}")
        return matrix_df
    
    def searchOrtholog(self, tax, orthologs_dataset, proteins):
        logging.debug(f"Searching orthologs for taxon {tax}...")
        tax_orthologs = orthologs_dataset[orthologs_dataset[f"{self._reference_species_column}_taxID"]==tax]
        if not tax_orthologs.empty:
            ort_counts = []
            orthologs_found = 0
            for protein in proteins:
                query = tax_orthologs[f"{self._reference_species_column}"].where(tax_orthologs[self._proteome_column]==protein).dropna() 
                if query.empty:
                    ort_counts.append(0)
                else:
                    ort_counts.append(query.count()) #Replace count() with tolist() if want to save the proteins IDs instead of just the count
                    orthologs_found+=1
            logging.debug(f"Found {orthologs_found} proteins with at least one ortholog in taxon {tax}.")
            return (tax, ort_counts)
        else:
            logging.warning(f"No ortholog has been found for taxon {tax} within this proteome.")


    

class TaxaMapping():

    @staticmethod
    def mapUniprot2Taxid_NCBI(uniprotIDs: set, prot_accession2taxid: os.path = "./data/prot.accession2QfOtaxid.gz") -> Dict:
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
                                names=["UniprotKB-AC", "GO", "NCBI-taxon"], #TODO: drop GO column
                                usecols=["UniprotKB-AC", "NCBI-taxon"],
                                dtype={"UniprotKB-AC":"string", "NCBI-taxon": "int32"})
        uniprot2taxid = uniprot2taxid[uniprot2taxid["UniprotKB-AC"].isin(uniprotIDs)] 
        # Check if all the uniprotIDs have a corresponding taxID
        if (uniprot2taxid["UniprotKB-AC"].nunique() != len(uniprotIDs)):
            logging.warning(f"The tax IDs for {len(uniprotIDs)-uniprot2taxid['UniprotKB-AC'].nunique()} uniprotKB accession numbers couldn't be found.")
            #logging.debug(set.difference(uniprotIDs, set(uniprot2taxid["UniprotKB-AC"].unique())))
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
