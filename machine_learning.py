# Basic imports
import logging
import sys
import pdb
import argparse
import os
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


from helpers.helper_functions import filterOutByFrequency

class ML():

    mlb = MultiLabelBinarizer()

    def __init__(self, phylogenetic_profile_matrix: pd.DataFrame):
        self.training_matrix_X = phylogenetic_profile_matrix.iloc[:, :-1] #values
        self.training_matrix_Y = pd.DataFrame(self.mlb.fit_transform(phylogenetic_profile_matrix[phylogenetic_profile_matrix.columns[-1]]), columns=self.mlb.classes_, index=phylogenetic_profile_matrix.index)
        self.models_available = self.create_baseline_classifiers()

    @property
    def labels(self):
        return self.training_matrix_Y
    
    @property
    def features(self):
        return self.training_matrix_X

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
        free_ncores = N_USABLE_CORES-cv if N_USABLE_CORES>cv else -2 # -2 means all cores available except 1
        models = {}
        # Inherently multilabel
        models['Dummy'] = DummyClassifier(random_state=seed, strategy='prior')
        models['RandomForest'] = RandomForestClassifier(random_state=seed, n_jobs=free_ncores)
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
        logging.debug(f"The following models are available: {', '.join(self.models_available.keys())}")
        logging.debug(f"Shape of X training matrix: {self.training_matrix_X.shape}")
        logging.debug(f"Shape of Y training matrix: {self.training_matrix_Y.shape}")
        logging.info("Assessing models...")
        summary = pd.DataFrame()
        for name, model in self.models_available.items():
            logging.info(f"Crossvalidating {name} model...")
            result = pd.DataFrame(cross_validate(model, self.training_matrix_X, self.training_matrix_Y, cv=cv, scoring=metrics, n_jobs=cv if N_USABLE_CORES>cv else 1))
            mean = result.mean().rename('{}_mean'.format)
            std = result.std().rename('{}_std'.format)
            summary[name] = pd.concat([mean, std], axis=0)
        return summary.sort_index()

    def train_model(self, model_name: str):
        model = self.models_available.get(model_name)
        if model: return model.fit()
        else:
            raise Exception("Specified model is not available")
    

def runML(pp_matrix: pd.DataFrame, min: int = None, max: int = None, results_file: os.path = ""):
    pp_matrix_training = pp_matrix[pp_matrix["GO_IDs"].str.len()>0] #the training dataset must be labeled
    if min: logging.info(f"Filtering out GO terms with less than {min} ocurrences.")
    if max: logging.info(f"Filtering out GO terms with more than {max} ocurrences.")
    pp_matrix_training = pp_matrix_training.assign(GO_IDs=filterOutByFrequency(pp_matrix_training["GO_IDs"], min_threshold=min, max_threshold=max))
    logging.info(f"Shrinked number of GO terms from {pp_matrix['GO_IDs'].explode().unique().size} to {pp_matrix_training['GO_IDs'].explode().unique().size}.")
    ml = ML(pp_matrix_training)
    assess_summary = ml.assess_models()
    logging.info(f"Assess summary:\n {assess_summary}")
    if results_file:
        assess_summary.to_csv(results_file, sep="\t", header=True, index=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Perform only Machine Learning", epilog="Enjoy!")
    parser.add_argument("--pp-matrix", type=str, required=False, default="drive/MyDrive/TFG/phylogenetic_profile_matrix_pres-absc_v2.tab", help="Phylogenetic Profiling matrix")
    parser.add_argument("--min-gos", type=int, required=False, default=None, help="Min number of GO terms' ocurrences,")
    parser.add_argument("--max-gos", type=int, required=False, default=None, help="Max number of GO terms' ocurrences,")
    parser.add_argument("--ml-results", type=str, required=False, default="drive/MyDrive/TFG/MTP_last-ML_assesment.tab", help="Filename for the Machine Learning models assessment results.")
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
                                header=0, index_col=0,  
                                converters={"GO_IDs": lambda x:  list(filter(None, x.split(",")))}) # if we don't filter there are no empty lists but lists with empty strings: [''] (its lenght is 1, not 0))
    runML(pp_matrix, min=args.min_gos, max=args.max_gos, results_file=args.ml_results)


