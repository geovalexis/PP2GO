# Basic imports
import logging
import sys
import pdb
import argparse
import os
from numpy.core.fromnumeric import mean, size
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
        self.matrix_Y = phylogenetic_profile_matrix.iloc[:, -1]
        self.one_hot_encoded_matrix_Y = pd.DataFrame(self.mlb.fit_transform(self.matrix_Y), columns=self.mlb.classes_, index=phylogenetic_profile_matrix.index)
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
    def create_baseline_classifiers(seed=8, cv=5):
        """Create a list of baseline classifiers.
        
        Parameters
        ----------
        seed: (optional) An integer to set seed for reproducibility
        Returns
        -------
        A list containing tuple of model's name and object.
        """
        free_ncores = N_USABLE_CORES-cv if N_USABLE_CORES>cv else -2    # -2 means all cores available except 1
        models = {}
        # Inherently multilabel
        models['Dummy'] = DummyClassifier(random_state=seed, strategy='prior')
        models['RandomForest'] = RandomForestClassifier(random_state=seed, n_jobs=1) # Random Forest consumes a lot of memory to be run in parallel
        models['KNN'] = KNeighborsClassifier(n_jobs=free_ncores)
        models['NeuralNetwork'] = MLPClassifier(random_state=seed)

        # No support for multilabel unless using OneVSRestClassifier or ClassifierChain
        models['SupportVectorMachine'] = OneVsRestClassifier(SVC(random_state=seed, probability=True), n_jobs=free_ncores)
        models['GradientBoosting'] = OneVsRestClassifier(GradientBoostingClassifier(random_state=seed), n_jobs=free_ncores)
        models['MultinomialNB'] = OneVsRestClassifier(MultinomialNB(), n_jobs=free_ncores)
        return models

    def assess_models(self, cv=5, metrics=("accuracy", "f1_macro")):
        """Provide summary of cross validation results for models.
        
        Parameters
        ----------
        X: A pandas DataFrame containing feature matrix
        y: A pandas Series containing target vector
        models: A list of models to train
        cv: (optional) An integer to set number of folds in cross-validation
        metrics: (optional) A list of scoring metrics or a string for a metric
        Returns
        -------
        A pandas DataFrame containing summary of baseline models' performance.
        """
        logger.debug(f"The following models are available: {', '.join(self.models_available.keys())}")
        logger.debug(f"Shape of X training matrix: {self.matrix_X.shape}")
        logger.debug(f"Shape of Y training matrix: {self.one_hot_encoded_matrix_Y.shape}")
        logger.info("Assessing models...")
        summary = pd.DataFrame()
        for name, model in self.models_available.items():
            logger.info(f"Crossvalidating {name} model...")
            result = pd.DataFrame(cross_validate(model, self.matrix_X, self.one_hot_encoded_matrix_Y, cv=cv, scoring=metrics, n_jobs=cv if N_USABLE_CORES>cv else 1))
            mean = result.mean().rename('{}_mean'.format)
            for m in metrics:
                logger.debug(f"Mean {m}: {mean['test_'+m+'_mean']}")
            std = result.std().rename('{}_std'.format)
            summary[name] = pd.concat([mean, std], axis=0)
        return summary.sort_index()

    def train_model(self, model_name: str):
        model = self.models_available.get(model_name)
        if model: return model.fit()
        else:
            raise Exception("Specified model is not available")

    def one_vs_rest(self, model_name: str, go_term: str, resampling: int):
        model = self.models_available.get(model_name)
        positive_samples = self.one_hot_encoded_matrix_Y[self.one_hot_encoded_matrix_Y[go_term] == 1].index
        negative_samples = self.one_hot_encoded_matrix_Y[self.one_hot_encoded_matrix_Y[go_term] == 0].index
        results = pd.DataFrame()
        logging.debug(f"Assessing {model_name} for GO term {go_term}...")
        for i in range(resampling):
            negative_samples_subset = np.random.choice(negative_samples, size=len(positive_samples)) #NOTE: positive and negative sample must be of the same size to avoid class imbalance
            training_samples = positive_samples.union(negative_samples_subset)
            X_train = self.matrix_X.loc[training_samples.tolist(), :]
            y_train = self.one_hot_encoded_matrix_Y.loc[training_samples.tolist(), go_term]
            results[i] = pd.DataFrame(cross_validate(model, X_train, y_train, cv=5, scoring=("accuracy", "f1_macro", "roc_auc"), n_jobs=-1)).mean()
        mean_results = results.mean(axis=1)
        logging.debug(f"Results after resampling {resampling} times:\n{mean_results}")
        return mean_results


def model_selection(pp_matrix: pd.DataFrame, min: int = None, max: int = None, results_file: os.path = ""):
    #pp_matrix_training = pp_matrix[pp_matrix["GO_IDs"].str.len()>0] #the training dataset must all be labeled
    pp_matrix_training = pp_matrix
    pp_matrix_training_size_before = pp_matrix_training['GO_IDs'].explode().unique().size

    if min: 
        logger.info(f"Filtering out GO terms with less than {min} ocurrences.")
        pp_matrix_training = pp_matrix_training.assign(GO_IDs=filterOutByMinFrequency(pp_matrix_training["GO_IDs"], min_threshold=min)).dropna() #NOTE: very import to drop those without any value

    if max: 
        logger.info(f"Filtering out GO terms with more than {max} ocurrences.")
        pp_matrix_training = pp_matrix_training.assign(GO_IDs=filterOutByMaxFrequency(pp_matrix_training["GO_IDs"], max_threshold=max)).dropna() #NOTE: very import to drop those without any value    

    logger.info(f"Filtering out GO terms present in all samples..")
    pp_matrix_training = pp_matrix_training.assign(GO_IDs=filterOutByExactFrequency(pp_matrix_training["GO_IDs"], freq=pp_matrix_training.shape[0])).dropna() # Drop those GO terms that are present in all samples (not informative and induces bias)
    logger.info(f"Shrinked number of distinct annotated GO terms from {pp_matrix_training_size_before} to {pp_matrix_training['GO_IDs'].explode().unique().size}.")
    labels_len_per_protein = pp_matrix_training["GO_IDs"].str.len()
    logger.info(f"Ocurrences mode: {labels_len_per_protein.mode()[0]}")
    logger.info(f"Ocurrences median: {labels_len_per_protein.median()}")
    logger.info(f"Ocurrences mean: {labels_len_per_protein.mean()}")
    #pp_matrix_training = pp_matrix_training[pp_matrix_training["GO_IDs"].str.len()>labels_len_mode]
    ml = ML(pp_matrix_training)
    assess_summary = ml.assess_models()
    logging.info(f"Assess summary:\n {assess_summary}")
    if results_file:
        assess_summary.to_csv(results_file, sep="\t", header=True, index=True)

def one_vs_rest_assessment(pp_matrix: pd.DataFrame, model_name: str, GO_terms: list, resampling_size: int = 1):
    pp_matrix_training = pp_matrix.loc[9606]
    ml = ML(pp_matrix_training)
    #logging.info(f"Assess summary:\n {assess_summary}")
    results = pd.DataFrame()
    for go_term in GO_terms:
        assess_summary = ml.one_vs_rest(model_name, go_term, resampling_size)
        results[go_term] = assess_summary
    return results

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Perform only Machine Learning", epilog="Enjoy!")
    parser.add_argument("--pp-matrix", type=str, required=False, default="results/MTP/MtP_201912-pp_matrix_counts.tab", help="Phylogenetic Profiling matrix")
    parser.add_argument("--min-gos", type=int, required=False, default=None, help="Min number of GO terms' ocurrences,")
    parser.add_argument("--max-gos", type=int, required=False, default=None, help="Max number of GO terms' ocurrences,")
    parser.add_argument("--ml-results", type=str, required=False, default="MTP_201912-prova-ML_assesment.tab", help="Filename for the Machine Learning models assessment results.")
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
    #model_selection(pp_matrix, min=args.min_gos, max=args.max_gos, results_file=args.ml_results)
    one_vs_rest_assessment(pp_matrix, "RandomForest", ["GO:0005739"], resampling_size=5)

