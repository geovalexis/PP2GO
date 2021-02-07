import pdb
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
        pdb.set_trace()
        return pd.get_dummies(column, 
                            prefix="", prefix_sep=""
                            ).reset_index().groupby("index").first()

if __name__ == "__main__":
    pp_matrix = pd.read_table("drive/MyDrive/TFG/phylogenetic_profile_matrix_pres-absc_v2.tab", header=0, index_col=0)
    pp_matrix_training = pp_matrix.dropna(subset=["GO_IDs"])
    print(ML(pp_matrix_training).getLabels)