# Phylogenetic Profiling to Gene Ontology (PP2GO)
Machine Learning to predict functional relationships of uncharacterized proteins through on orthology and co-evolution information.
This project has been developed as a part of my Final Bachelor Thesis for the double degree in Biotechnology and Computer Engineering at URV (Tarragona, Spain).
Both methodology and resources have been provided by the Comparative Genomics lead by Toni Gabaldon.
For more information about the methodology of this pipeline, please refer to the [wiki](https://github.com/geovalexis/PP2GO/wiki) section of the repo.

## Repository tree
  *	**data/**: various data needed for the pipeline to run. They are not enough to run the pipeline. GitHub strictly limits the size of individual files to 100MB, so some files such as the GOA annotations or taxid_mapping could not be uploaded. These must be downloaded manually. See the instructions on 1_Data_Retrieval.ipynb.
  *	**orthologs/**: folder with the orthologs datasets from different sources. As before, the size limitation does not allow us to upload the entire MetaPhOrs’ orthologs dataset but only a subset of it, which can be used for testing purposes. 
  *	**notebooks/**: jupyter notebooks (powered by Google Colab) used as a testing/debugging environment for the algorithms. It also includes the implementation of all the plots shown in the work as well as most of the data analysis performed to achieve the results. 
  *	**src/**: source code of the pipeline used to test different filters and approaches. It was the code optimized to be run in the MareNostrum4. 
  *	**pp2go_env.yml**: pipeline’s package requirements to install with Anaconda. 

## 1. Installation
For being able to install this pipeline, you will need to have [Anaconda](https://docs.anaconda.com/anaconda/install/) installed . It is highly recommended that you create a new environment before installing the dependencies. You can do it directly by importing the environment provided:
```bash
conda env create -f pp2go_env.yml
conda activate pp2go
```
## 2. Datasets
The main script of the pipeline is found at **src/pp2go.py**. For correctly run the pipeline, you will need the following files:
* IdMapping File, from the official [Uniprot FTP server](http://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/idmapping/). 
* GAF file, from the official [EBI FTP sever](http://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/)
* OBO file of Gene Ontology, but this is already found at **data/go.obo**. if not provided, it will be automatically downloaded from the the official [GO consortium server](http://geneontology.org/docs/download-ontology/) and moved to **data/** folder.
* Orthologs dataset (a mock dataset is provided at **orthologs/MtP**)


## 3. Running pipeline
There are different parameters that can be fill, but a normal run would be as follows:
```bash
python src/pp2go.py \
  --idmapping_file ./data/idmapping_selected_qfo_subset.tab.gz \
  --gaf-file ./data/goa_uniprot_qfo.gaf.gz \
  --orthologs ./orthologs/MTP/MtP_201912.tab  \
  --include-go-parents \
  --min-gos 100 \
  --max-gos 5000 \
  --pp-matrix ./results/MTP/MtP_201912-include_go_parents-min_gos_100-max_gos_5000-pp_matrix.tab  \
  --ml-results ./results/MTP/MtP_201912-include_go_parents-min_gos_100-max_gos_5000-ML_assessment.tab  \
```
## 4. Output files description
* PP or Phylogenetic Profiling matrix: it consists of a tabular with the phylogenetic profiling matrix when all GO terms have been assigned and with all the specified filters applied (minimum/maximum number of GO terms, evidence codes, etc).
* ML or Machine Learning assessment: it consists of a tabular file with the performance scores of the Machine Learning analysis (either the allVSall or oneVSrest). 

## Parameters complete list
```
  -h, --help            show this help message and exit
  --gaf-file GAF_FILE   Gene Ontology annotation file.
  --idmapping_file IDMAPPING_FILE
                        Identifiers mapping file from Uniprot.
  --orthologs ORTHOLOGS
                        Orthologs input.
  --filter-by-sp        Select if wants to filter by only Swiss Prot proteins
  --pres-abs            Compute presence/abscense PP Matrix instead of counts
                        PP matrix.
  --proteome-species [PROTEOME_SPECIES [PROTEOME_SPECIES ...]]
                        Space separated list of species whose proteins will be
                        used for the Phylogenetic Profiling Matrix. Human
                        proteome will be taken by default.
  --proteome PROTEOME   File that contains a line-separated list of proteins
                        to be take as input proteome.
  --reference-species [REFERENCE_SPECIES [REFERENCE_SPECIES ...]]
                        Space separated list of reference organisms on which
                        the orthologs will be searched for. By default all
                        available will be taken.
  --pp-matrix PP_MATRIX
                        Name of the Phylogenetic Profiling Matrix if wants to
                        be saved.
  --go-aspects [{P,C,F} [{P,C,F} ...]]
                        GO aspect/ontology. By default only Biological Process
                        will be taken.
  --set-as-root SET_AS_ROOT
                        Set a given GO term as root, so only their children
                        can be assigned
  --include-go-parents  Include all the lineage for each GO term assign to a
                        protein.
  --min-level MIN_LEVEL
                        Minimum level of GO terms that can be assigned
  --max-level MAX_LEVEL
                        Maximum level of GO terms that can be assigned
  --min-gos MIN_GOS     Min number of GO terms' ocurrences,
  --max-gos MAX_GOS     Max number of GO terms' ocurrences,
  --ml-mode {all,one}   Mode of the ML algorithm. Choices: multilabel (all) or
                        binary classification (one).
  --ml-results ML_RESULTS
                        Filename for the Machine Learning models assessment
                        results.
  -v, --verbose         Verbose logging.
```

