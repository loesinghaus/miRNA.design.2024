import pandas as pd
from .mirna_combinations import *

def additive_expression(expression: pd.DataFrame, 
                        combinations: list) -> pd.DataFrame:
    """This function takes a dataframe of gene expression and a list of miRNA combinations
    and returns a dataframe with the expression of the miRNAs in the combinations
    using an additive model."""
    
    # create a new dataframe
    expression_comb = {}
    # invert the expression values
    expression = 1/expression
    for comb in combinations:
        expression_comb[comb] = (expression.loc[comb,:]-1).sum(axis=0)+1
    # invert the combined expression values
    expression_comb = 1/pd.DataFrame(expression_comb)
    return expression_comb.T

def add_knockdown(knockdown, combinations, design_miRNA_dict):
    """This functions takes knockdown (expression) values caused by a single microRNA, microRNA combinations, and a dictionary
    that maps microRNA combinations to design names. It returns the expression values caused by the microRNA combinations.
    
    knockdown: dataframe with knockdown/expression values caused by a single microRNA (index: microRNA names, columns: cell lines)
    combinations: list of tuples with microRNA names
    design_miRNA_dict: dictionary that maps microRNA combinations to design names."""

    add = additive_expression(knockdown, combinations)

    new_index = [design_miRNA_dict[tuple(sorted(mirnas))] for mirnas in add.index.to_list()]
    add.index = new_index

    return add

def add_mirna_expression(mirna_expr, construct_df):
    """Takes a dataframe with microRNA expressions and a dataframe with construct information.
    Returns a dataframe with the added expression of the microRNAs in the constructs.
    
    mirna_expr: dataframe with microRNA expression. This is assumed to be linear. (index: microRNA names, columns: cell lines)
    construct_df: dataframe with construct information (index: construct names, columns: microRNA names)"""

    combs = get_combinations(construct_df)
    added_expression = pd.DataFrame(columns=mirna_expr.columns, index=construct_df.index)

    for i, comb in enumerate(combs):
        added_expression.iloc[i,:] = mirna_expr.loc[comb,:].sum(axis=0)

    return added_expression

def max_mirna_expression(mirna_expr, construct_df):
    """Takes a dataframe with microRNA expressions and a dataframe with construct information.
    Returns a dataframe with the added expression of the microRNAs in the constructs.
    
    mirna_expr: dataframe with microRNA expression (index: microRNA names, columns: cell lines)
    construct_df: dataframe with construct information (index: construct names, columns: microRNA names)"""

    combs = get_combinations(construct_df)
    added_expression = pd.DataFrame(columns=mirna_expr.columns, index=construct_df.index)

    for i, comb in enumerate(combs):
        added_expression.iloc[i,:] = mirna_expr.loc[comb,:].max(axis=0)

    return added_expression

def add_mirna_combs(mirna_expr, combs):
    """Takes a dataframe with microRNA expressions and tuples of combinations.
    Returns a dataframe with the added expression of the microRNAs in the constructs.
    
    mirna_expr: dataframe with microRNA expression (index: microRNA names, columns: cell lines)
    construct_df: dataframe with construct information (index: construct names, columns: microRNA names)"""
    # if isinstance(combs[0], tuple):
    multiindex = pd.MultiIndex.from_tuples(combs, names=[f'miRNA{i+1}' for i in range(len(combs[0]))])
    added_expression = pd.DataFrame(columns=mirna_expr.columns, index=multiindex)
    added_expression = added_expression.astype("float")
    # else:
    #     added_expression = pd.DataFrame(columns=mirna_expr.columns, index=combs)

    for i, comb in enumerate(combs):
        added_expression.loc[comb,:] = mirna_expr.loc[comb,:].sum(axis=0).values

    return added_expression

def max_mirna_combs(mirna_expr, combs):
    """Takes a dataframe with microRNA expressions and tuples of combinations.
    Returns a dataframe with the added expression of the microRNAs in the constructs.
    
    mirna_expr: dataframe with microRNA expression (index: microRNA names, columns: cell lines)
    construct_df: dataframe with construct information (index: construct names, columns: microRNA names)"""
    # if isinstance(combs[0], tuple):
    multiindex = pd.MultiIndex.from_tuples(combs, names=[f'miRNA{i}' for i in range(len(combs[0]))])
    added_expression = pd.DataFrame(columns=mirna_expr.columns, index=multiindex)
    added_expression = added_expression.astype("float")
    # else:
    #     added_expression = pd.DataFrame(columns=mirna_expr.columns, index=combs)

    for i, comb in enumerate(combs):
        added_expression.loc[comb,:] = mirna_expr.loc[comb,:].max(axis=0).values

    return added_expression