import pdb
from numpy.lib.arraysetops import setdiff1d
import pandas as pd
import numpy as np


class ML():

    def __init__(self, phylogenetic_profile_matrix: pd.DataFrame):
        self.training_matrix_X = phylogenetic_profile_matrix.iloc[:, :-1] #values
        self.training_matrix_Y = self.applyOneHotEncoding(phylogenetic_profile_matrix[phylogenetic_profile_matrix.columns[-1]]) #labels

    @property
    def getLabels(self):
        return self.training_matrix_Y
    
    @property
    def getValues(self):
        return self.training_matrix_X

    @staticmethod
    def applyOneHotEncoding(column: pd.Series):
        return pd.get_dummies(column.explode(), 
                            prefix="", prefix_sep=""
                            ).reset_index().groupby("index").first()

def filterOutByFrequency(GO_ids: pd.Series, min_threshold: int = None, max_threshold: int = None):
    go_ids_counts = GO_ids.explode().value_counts()
    min_threshold = min_threshold if min_threshold else go_ids_counts.max()
    max_threshold = max_threshold if max_threshold else go_ids_counts.min()
    out_values = go_ids_counts[(go_ids_counts<min_threshold) | (go_ids_counts>max_threshold)].index
    GO_ids_filtered = GO_ids.apply(lambda x, y=out_values: np.setdiff1d(x, y)).copy()
    return GO_ids_filtered
    

if __name__ == "__main__":
    pp_matrix = pd.read_table("drive/MyDrive/TFG/phylogenetic_profile_matrix_pres-absc_v2.tab", 
                                header=0, index_col=0,  
                                converters={"GO_IDs": lambda x:  list(filter(None, x.split(",")))}) # if we don't filter there are no empty lists but lists with empty strings: [''] (its lenght is 1, not 0))
    pp_matrix_training = pp_matrix[pp_matrix["GO_IDs"].str.len()>0]
    pp_matrix_training = pp_matrix_training.assign(GO_IDs=filterOutByFrequency(pp_matrix_training["GO_IDs"], min_threshold=100, max_threshold=1000))
    print(ML(pp_matrix_training).getLabels)