import os
import logging
from typing import Dict
import pandas as pd
import numpy as np
import pdb
import sys
from joblib import Parallel, delayed

__author__ = "Geovanny Risco"
__email__ = "geovanny.risco@bsc.es"
__version__ = "0.1"

logger = logging.getLogger(__name__)


class PhylogeneticProfiling():
    """ 
    It takes an orthologs dataset as input in the following format: uniprotid1 | uniprotid2.
    Then it maps the corresponding tax Id for each protein (see TaxaMapping class), obtaining a matrix with:  uniprotid1 | uniprotid2 | uniprotid1_taxID | uniprot2_taxID. 
    Both set of proteins are supposed to be interexchangable, so in this case the first set of uniprot IDs will act as proteome 
    dataset and the second as orthologs dataset. 
    """

    def __init__(self, idmapping_file: os.path, orthologs: pd.DataFrame, reference_species: list = [], onProteins: list = [], onSpecies: list = []):
        self._proteome_column = orthologs.columns[0]
        self._reference_species_column = orthologs.columns[1]
        self._species_availables = None
        self._orthologs = self.mapAndfilterOrthologs(idmapping_file, orthologs, reference_species, onProteins, onSpecies)

    def mapAndfilterOrthologs(self, idmapping_file, orthologs, reference_species, onProteins, onSpecies):
        """ Map taxIDs to the orthologs and filter dataset by the given reference species and proteins/species
        """
        taxd_orthologs = self.mapTaxIDs(idmapping_file, orthologs, onColumns=[self._proteome_column, self._reference_species_column], dropUnmatched=True) #TODO: Support for OrthoXML format
        self._species_availables = list(taxd_orthologs[f"{self._reference_species_column}_taxID"].unique()) #Save all species available in the dataset
        if reference_species:
            filtered_orthologs_dataset = self.filterByReferenceSpecies(taxd_orthologs, set(reference_species))      
        else:
            filtered_orthologs_dataset = taxd_orthologs  

        if onProteins:
            filtered_orthologs_dataset = self.filterByProteins(taxd_orthologs, set(onProteins))
        elif onSpecies:
            filtered_orthologs_dataset = self.filterBySpecies(taxd_orthologs,  set(onSpecies))
        return filtered_orthologs_dataset

    def filterByReferenceSpecies(self, taxd_orthologs_dataset, reference_species):
        logger.debug("Filtering orthologs dataset by reference species...")
        return taxd_orthologs_dataset[taxd_orthologs_dataset[f"{self._reference_species_column}_taxID"].isin(reference_species)]

    def filterByProteins(self, taxd_orthologs_dataset, onProteins):
        logger.debug("Filtering orthologs dataset by proteins...")
        return taxd_orthologs_dataset[taxd_orthologs_dataset[self._proteome_column].isin(onProteins)]

    def filterBySpecies(self, taxd_orthologs_dataset, onSpecies):
        logger.debug("Filtering orthologs dataset by species...")
        return taxd_orthologs_dataset[taxd_orthologs_dataset[f"{self._proteome_column}_taxID"].isin(onSpecies)]

    @property
    def referenceSpecies(self):
        return self._species_availables
    
    @property
    def onSpecies(self):
        return list(self._orthologs[f"{self._proteome_column}_taxID"].unique())

    @property
    def onProteins(self):
        return list(self._orthologs[self._proteome_column].unique())
    

    @staticmethod
    def mapTaxIDs(idmapping: os.path, df: pd.DataFrame, onColumns: list, dropUnmatched: bool = True):
        """ Maps the taxIds for each uniprotId on the given columns of the dataframe.

        Args:
            df (pd.DataFrame): dataframe to map
            onColumns (List[str]): columns of the dataframe where the uniprotids are
            dropUnmatched (bool): whether drop o not those rows that did not have a match (no taxon assigned)
        Returns:
            df: same dataframe as passed in but with the taxIDs columns added with 
                the following format <onColumn>_taxID
        """
        all_uniprotids = []
        for column in onColumns:
            all_uniprotids.extend(df[column].unique())
        logger.info("Starting taxID assignment...")
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
        logger.debug(f"Total number of unique proteins: {len(proteome)}")
        logger.info("Computing Phylogenetic Profiling matrix...")
        with Parallel(n_jobs=-1) as run_in_parallel:
            matrix = dict(filter(None, 
                        run_in_parallel(delayed(
                            self.searchOrtholog) (tax, self._orthologs, proteome, self._proteome_column, self._reference_species_column) 
                                for tax in self.referenceSpecies)))
        logger.info(f"...found orthologs in {len(matrix)} out of {len(self._species_availables)} taxons.")
        matrix_df = pd.DataFrame(matrix, index=[proteome_taxa,proteome])
        logger.info(f"Final shape of the matrix: {matrix_df.shape}")
        return matrix_df
    
    @staticmethod
    def searchOrtholog(tax, orthologs_dataset, proteins, proteome_column, reference_species_column):
        logger.debug(f"Searching orthologs for taxon {tax}...")
        tax_orthologs = orthologs_dataset[orthologs_dataset[f"{reference_species_column}_taxID"]==tax]
        if not tax_orthologs.empty:
            ort_counts = []
            orthologs_found = 0
            for protein in proteins:
                query = tax_orthologs[f"{reference_species_column}"].where(tax_orthologs[proteome_column]==protein).dropna() 
                if query.empty:
                    ort_counts.append(0)
                else:
                    ort_counts.append(query.count()) #Replace count() with tolist() if want to save the proteins IDs instead of just the count
                    orthologs_found+=1
            logger.debug(f"Found {orthologs_found} proteins with at least one ortholog in taxon {tax}.")
            return (tax, ort_counts)
        else:
            logger.warning(f"No ortholog has been found for taxon {tax} within this proteome.")


    

class TaxaMapping():

    @staticmethod
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
            logger.warning(f"The tax ID for {len(uniprotIDs)-uniprot2taxid['accession'].nunique()} uniprotKB accession numbers couldn't be found.")
            logger.debug(set.difference(uniprotIDs, set(uniprot2taxid["accession"].unique())))
        return uniprot2taxid.set_index("accession")["taxid"].to_dict()

    
    @staticmethod
    def mapUniprot2Taxid_Uniprot(uniprotIDs: set, idmapping: os.path) -> Dict: 
        """ Translate a set of genes to its corresponding taxIDs given their uniprotKB accession numbers.

        IMPORTANT:  This code will not work with the original `idmapping_selected_subset.tab.gz` file.
                    See '1_Data_Retrieval' notebook, section 'Id mapping file', to learn how this idmapping 
                    file was downloaded and then filtered in order to save memory and disk storage.
        Args:
            uniprotIDs (set): set of genes to translate
            idmapping (os.path): path to the idmapping_selected_filtered.tab gzip compressed file 
                                from the official Uniprot ftp server.
                                See ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/idmapping/ 
        """    
        uniprot2taxid = pd.read_table(idmapping,
                                compression="gzip", 
                                names=["UniprotKB-AC", "GO", "NCBI-taxon"], 
                                usecols=["UniprotKB-AC", "NCBI-taxon"], #NOTE: GO is assigned with the GAF file as it allows more filtering
                                dtype={"UniprotKB-AC":"string", "NCBI-taxon": "int32"})
        uniprot2taxid = uniprot2taxid[uniprot2taxid["UniprotKB-AC"].isin(uniprotIDs)] 
        # Check if all the uniprotIDs have a corresponding taxID
        if (uniprot2taxid["UniprotKB-AC"].nunique() != len(uniprotIDs)):
            logger.warning(f"The tax IDs for {len(uniprotIDs)-uniprot2taxid['UniprotKB-AC'].nunique()} uniprotKB accession numbers couldn't be found.")
            #logger.debug(set.difference(uniprotIDs, set(uniprot2taxid["UniprotKB-AC"].unique()))) #NOTE: if the differences are too high the whole text can be annoying 
        return uniprot2taxid.set_index("UniprotKB-AC")["NCBI-taxon"].to_dict()

                