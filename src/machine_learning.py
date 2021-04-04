# Basic imports
import logging
import sys
import pdb
import argparse
import pandas as pd
import numpy as np

#For multiprocessing
from os import sched_getaffinity
N_USABLE_CORES = len(sched_getaffinity(0)) # This give us the actual number of usable cores, unlike cpu_count(), 
                                           # which give us the total (in MN4 would always say 48, although we had assigned a cpus-per-task of 24)

# Scikit learn utils 
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.multiclass import OneVsRestClassifier
from sklearn.multioutput import ClassifierChain
from sklearn.metrics import accuracy_score
from sklearn.model_selection import cross_validate

# Models that support multilabel by default
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier

# Models that do not inherently support multilabel (A layer like OneVsRestClassifier must be used)
from sklearn.svm import SVC
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.naive_bayes import MultinomialNB

from helpers.helper_functions import filterOutByExactFrequency, filterOutByMaxFrequency, filterOutByMinFrequency

__author__ = "Geovanny Risco"
__email__ = "geovanny.risco@bsc.es"
__version__ = "0.1"

logger = logging.getLogger(__name__)

class ML():

    mlb = MultiLabelBinarizer()

    def __init__(self, phylogenetic_profile_matrix: pd.DataFrame):
        self.matrix_X = phylogenetic_profile_matrix.iloc[:, :-1] 
        self.matrix_Y = pd.DataFrame(self.mlb.fit_transform(phylogenetic_profile_matrix.iloc[:, -1]), columns=self.mlb.classes_, index=phylogenetic_profile_matrix.index)
        self.models_available = self.create_baseline_classifiers()

    @property
    def labels(self):
        return self.matrix_Y
    
    @property
    def features(self):
        return self.matrix_X

    @property
    def modelsAvailable(self):
        return self.models_available

    # From https://towardsdatascience.com/simple-way-to-find-a-suitable-algorithm-for-your-data-in-scikit-learn-python-9a9710c7c0fe
    @staticmethod
    def create_baseline_classifiers(seed=8, n_jobs=1):
        """Create a list of baseline classifiers.
        
        Args:
            seed (int, optional): an integer to set seed for reproducibility. Defaults to 8.
            cv (int, optional): k-fold of the cross-validation. Needed to calculate the number of free cores to use in the parallelization.
        Returns:
            models (list): a list containing tuple of model's name and object.
        """
        models = {}
        # Inherently multilabel
        models['Dummy'] = DummyClassifier(random_state=seed, strategy='prior')
        models['RandomForest'] = RandomForestClassifier(random_state=seed, n_jobs=1) # Random Forest consumes a lot of memory to be run in parallel, so better keep it with 1 core always
        models['KNN'] = KNeighborsClassifier(n_jobs=n_jobs)
        models['NeuralNetwork'] = MLPClassifier(random_state=seed)

        # No support for multilabel unless using OneVSRestClassifier or ClassifierChain
        models['SupportVectorMachine'] = OneVsRestClassifier(SVC(random_state=seed, probability=True), n_jobs=n_jobs)
        models['GradientBoosting'] = OneVsRestClassifier(GradientBoostingClassifier(random_state=seed), n_jobs=n_jobs)
        models['MultinomialNB'] = OneVsRestClassifier(MultinomialNB(), n_jobs=n_jobs)
        return models

    def assess_all_labels(self, models: list, cv=5, metrics=("accuracy", "f1_macro")):
        """ Assess all the labels together as a multilabel problem and provides a summary of the 
        cross validation results for the requested models.

        Args:
            models (list): A list of models to train
            cv (int, optional): An integer to set number of folds in cross-validation. Defaults to 5.
            metrics (tuple, optional): A list of scoring metrics or a string for a metric. Defaults to ("accuracy", "f1_macro").

        Returns:
            summary (pd.DataFrame): a pandas DataFrame containing a summary of the models' performance in the requested metrics.
        """
        models_to_assess = {model_name:self.models_available[model_name] for model_name in models}
        logger.debug(f"The following models will be evaluated: {', '.join(models_to_assess)}")
        logger.debug(f"Shape of X training matrix: {self.matrix_X.shape}")
        logger.debug(f"Shape of Y training matrix: {self.matrix_Y.shape}")
        logger.info("Assessing models...")
        summary = pd.DataFrame()
        for name, model in models_to_assess.items():
            logger.info(f"Crossvalidating {name} model...")
            result = pd.DataFrame(cross_validate(model, self.matrix_X, self.matrix_Y, cv=cv, scoring=metrics, n_jobs=-1))
            mean = result.mean().rename('{}_mean'.format)
            for m in metrics:
                logger.debug(f"Mean {m}: {mean['test_'+m+'_mean']}")
            std = result.std().rename('{}_std'.format)
            summary[name] = pd.concat([mean, std], axis=0)
        return summary.sort_index()

    def assess_one_label(self, go_term: str, resampling: int, models: list, cv=5, metrics=("accuracy", "f1_macro", "roc_auc")):
        """Assess only one label (binary classification problem), taking as many negative samples as positive sample (prevents class imbalance).

        Args:
            go_term (str): label name.
            resampling (int): number of times to resample (shuffle) the negative samples.
            models (list): names of model to assess
            cv (int, optional): K fold for Cross-Validation. Defaults to 5.
            metrics (tuple, optional): evaluation metrics to return. Defaults to ("accuracy", "f1_macro", "roc_auc").

        Returns:
            summary (pd.DataFrame): a pandas DataFrame containing a summary of the models' performance in the requested metrics.
        """
        models_to_assess = {model_name:self.models_available[model_name] for model_name in models}
        positive_samples = self.matrix_Y[self.matrix_Y[go_term] == 1].index
        negative_samples = self.matrix_Y[self.matrix_Y[go_term] == 0].index
        summary = pd.DataFrame()
        logging.info(f"Assessing GO term {go_term}...")
        for name, model in models_to_assess.items():
            logging.info(f"With {name}...")
            results = pd.DataFrame()
            for i in range(resampling):
                negative_samples_subset = np.random.choice(negative_samples, size=len(positive_samples)) #NOTE: positive and negative sample must be of the same size to avoid class imbalance
                training_samples = positive_samples.union(negative_samples_subset)
                X_train = self.matrix_X.loc[training_samples.tolist(), :]
                y_train = self.matrix_Y.loc[training_samples.tolist(), go_term]
                results[i] = pd.DataFrame(cross_validate(model, X_train, y_train, cv=cv, scoring=metrics, n_jobs=-1)).mean()
            results_mean = results.mean(axis=1).rename('{}_mean'.format)
            results_std = results.std(axis=1).rename('{}_std'.format)
            summary[name] = pd.concat([results_mean, results_std], axis=0)
            logging.debug(f"...results after resampling {resampling} times:\n{summary[name]}")
        return summary.sort_index()

def prepare_training_matrix(pp_matrix: pd.DataFrame, min: int, max: int):
    #pp_matrix = pp_matrix[pp_matrix["GO_IDs"].str.len()>0] #the training dataset must all be labeled
    pp_matrix_training = pp_matrix.copy()
    pp_matrix_training_size_before = pp_matrix['GO_IDs'].explode().unique().size

    if min: 
        logging.info(f"Filtering out GO terms with less than {min} ocurrences.")
        pp_matrix_training = pp_matrix_training.assign(GO_IDs=filterOutByMinFrequency(pp_matrix_training["GO_IDs"], min_threshold=min)).dropna() #NOTE: very import to drop those without any value

    if max: 
        logging.info(f"Filtering out GO terms with more than {max} ocurrences.")
        pp_matrix_training = pp_matrix_training.assign(GO_IDs=filterOutByMaxFrequency(pp_matrix_training["GO_IDs"], max_threshold=max)).dropna() #NOTE: very import to drop those without any value    

    logging.info(f"Filtering out GO terms present in all samples..")
    pp_matrix_training = pp_matrix_training.assign(GO_IDs=filterOutByExactFrequency(pp_matrix_training["GO_IDs"], freq=pp_matrix_training.shape[0])).dropna() # Drop those GO terms that are present in all samples (not informative and induces bias)
    logging.info(f"Shrinked number of distinct annotated GO terms from {pp_matrix_training_size_before} to {pp_matrix_training['GO_IDs'].explode().unique().size}.")
    labels_len_per_protein = pp_matrix_training["GO_IDs"].str.len()
    logging.info(f"Ocurrences mode: {labels_len_per_protein.mode()[0]}")
    logging.info(f"Ocurrences median: {labels_len_per_protein.median()}")
    logging.info(f"Ocurrences mean: {labels_len_per_protein.mean()}")
    #pp_matrix_training = pp_matrix_training[pp_matrix_training["GO_IDs"].str.len()>labels_len_mode]   
    return pp_matrix_training

def all_vs_all_assessment(pp_matrix_training: pd.DataFrame):
    ml = ML(pp_matrix_training)
    assess_summary = ml.assess_all_labels(models=ml.models_available.keys())
    logging.info(f"Assess summary:\n {assess_summary}")
    return assess_summary

def one_vs_rest_assessment(pp_matrix_training: pd.DataFrame, GO_terms: list, resampling_size: int = 1):
    ml = ML(pp_matrix_training)
    results = {}
    for go_term in GO_terms:
        assess_summary = ml.assess_one_label(go_term, resampling_size, models=["RandomForest"])
        results[go_term] = assess_summary
    results_df = pd.concat(results, axis=1)
    return results_df

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Perform only Machine Learning", epilog="Enjoy!")
    parser.add_argument("--pp-matrix", type=str, required=False, default="results/MTP/MtP_201912-pp_matrix_counts.tab", help="Phylogenetic Profiling matrix")
    parser.add_argument("--min-gos", type=int, required=False, default=None, help="Min number of GO terms' ocurrences,")
    parser.add_argument("--max-gos", type=int, required=False, default=None, help="Max number of GO terms' ocurrences,")
    parser.add_argument("--go-terms", type=argparse.FileType("rt"), required=False, default=None, help="File with a line-separated list of GO terms to assess individually.")
    parser.add_argument("--ml-results", type=str, required=False, default=None, help="Filename for the Machine Learning models assessment results.")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] -- %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stdout
    ) 
    pp_matrix = pd.read_table(args.pp_matrix, 
                                header=0, index_col=[0,1],  
                                converters={"GO_IDs": lambda x:  list(filter(None, x.split(",")))}) # if we don't filter there are no empty lists but lists with empty strings: [''] (its lenght is 1, not 0))
    pp_matrix_training = prepare_training_matrix(pp_matrix, min=args.min_gos, max=args.max_gos)
    results = one_vs_rest_assessment(pp_matrix_training, GO_terms=args.go_terms.read().splitlines() if args.go_terms else pp_matrix_training["GO_IDs"].explode().unique().tolist(), resampling_size=5)
    #results = all_vs_all_assessment(pp_matrix_training)
    if args.ml_results:
        results.to_csv(args.ml_results, sep="\t", header=True, index=True)

