import pdb
import pandas as pd
import numpy as np

def filterOutByFrequency(column: pd.Series, min_threshold: int = None, max_threshold: int = None) -> pd.Series:
    """ Take pandas' Series object and filter their values by a minimum and a maximun threshold.

    Args:
        column (pd.Series): Each row within the Series must be a list or any list-like structure (set, tuple, np.darray, etc.)
        min_threshold (int, optional): minumum number of times that a same value can be repeated. If not defined, it won't be applied.
        max_threshold (int, optional): maximum number of times that a same value can be repeated. If not defined, it won't be applied.

    Returns:
        column_filtered (pd.Series): copy of the Series with all the elements beyond the thresholds filtered out.
    """
    elements_counts = column.explode().value_counts()
    print(f"The element with the maximum number of ocurrences is {elements_counts.idxmax()} with {elements_counts.max()} ocurrences.", flush=True)
    print(f"The element with the maximum number of ocurrences is {elements_counts.idxmin()} with {elements_counts.min()} ocurrences.", flush=True)
    min_threshold = min_threshold if min_threshold else elements_counts.max()
    max_threshold = max_threshold if max_threshold else elements_counts.min()
    out_elements = elements_counts[(elements_counts<min_threshold) | (elements_counts>max_threshold)].index
    column_filtered = column.apply(lambda x, y=out_elements: np.setdiff1d(x, y)).copy()
    return column_filtered
