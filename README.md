# Phylogenetic Profiling to Gene Ontology (PP2GO)
Machine Learning to predict functional relationships of uncharacterized proteins through on orthology and co-evolution information.
This project has been developed as a Final Thesis for the double degree in Biotechnology and Computer Engineering at URV (Tarragona, Spain).
Both methodology and resources have been provided by the Comparative Genomics lead by Toni Gabaldon.

## Repository tree
  *	**data/**: various data needed for the pipeline to run. They are not enough to run the pipeline. GitHub strictly limits the size of individual files to 100MB, so some files such as the GOA annotations or taxid_mapping could not be uploaded. These must be downloaded manually. See the instructions on 1_Data_Retrieval.ipynb.
  *	**orthologs/**: folder with the orthologs datasets from different sources. As before, the size limitation does not allow us to upload the entire MetaPhOrs’ orthologs dataset but only a subset of it, which can be used for testing purposes. 
  *	**notebooks/**: jupyter notebooks (powered by Google Colab) used as a testing/debugging environment for the algorithms. It also includes the implementation of all the plots shown in the work as well as most of the data analysis performed to achieve the results. 
  *	**src/**: source code of the pipeline used to test different filters and approaches. It was the code optimized to be run in the MareNostrum4. 
  *	**pp2go_env.yml**: pipeline’s package requirements to install with Anaconda. 



