import numpy as np
import pandas as pd

def merge_identical_mirnas(df_expression, mirbase):
    """Merges identical miRNAs in the expression dataframe.
    Assumes log10 expression levels."""
    
    # First, we merge identical miRNAs
    mirbase["sequence_norm_short"] = mirbase["sequence_norm"].str[:18]

    # find out which miRNAs are identical
    mirna_groups = mirbase.groupby("sequence_norm_short").groups
        
    # add their expression levels
    for key, group in mirna_groups.items():
        # filter the group to those that are in the expression dataframe
        group = [x for x in group if x in df_expression.index]
        if len(group) > 1:
            # get the expression levels
            expression_levels = df_expression.loc[group, :]
            # make them linear
            expression_levels = 10**expression_levels
            # sum them
            expression_levels = np.log10(expression_levels.sum(axis=0))
            # set all mirnas in the group to the summed expression levels
            for mirna in group:
                for column in df_expression.columns:
                    # is the original value NaN? skip!
                    if np.isnan(df_expression.loc[mirna, column]):
                        continue
                    df_expression.loc[mirna, column] = expression_levels[column]
                
    return df_expression, mirna_groups