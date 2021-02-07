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
        return pd.get_dummies(column.explode(), 
                            prefix="", prefix_sep=""
                            ).reset_index().groupby("index").first()

if __name__ == "__main__":
    pp_matrix = pd.read_table("drive/MyDrive/TFG/phylogenetic_profile_matrix_pres-absc_v2.tab", 
                                header=0, index_col=0,  
                                converters={"GO_IDs": lambda x:  list(filter(None, x.split(",")))}) # if we don't filter there are no empty lists but lists with empty strings: [''] (its lenght is 1, not 0))
    pp_matrix_training = pp_matrix[pp_matrix["GO_IDs"].str.len()>0]
    print(ML(pp_matrix_training).getLabels)