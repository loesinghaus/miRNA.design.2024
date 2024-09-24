import pandas as pd
import numpy as np

# def levels_to_expression_elu_noscale(x, et_1=-0.15, et_2=-0.8, em_1=3.5, em_2=5, a3=3):
#     """This function assumes that miRNA repression is elu-like,
#     with an inflection at em_1 with an expression of et_1,
#     and a expression of et_2 at em_2. .
#     x is a numpy array of log10 expression values.
#     returns log10 expression values."""

#     # linear part
#     a1 = (et_2-et_1)/(em_2-em_1)
#     b1 = et_1-a1*em_1
#     linear = a1*x+b1

#     # exponential part
#     a2 = a1/a3
#     expo = a2*(np.exp(a3*(x-em_1))-1)+et_1

#     # adjust slightly
#     expo0 = a2*(np.exp(a3*(-em_1))-1)+et_1
#     linear = linear-expo0
#     expo = expo-expo0

#     return np.where(x<em_1, expo, linear)

# def levels_to_expression_elu_scales(datasets, et_1=-0.15, et_2=-0.8, em_1=3.5, em_2=5, a3=3, *scales):
#     # datasets is now a list of numpy arrays, scales is a list of scale values
#     results = []

#     for x, scale in zip(datasets, scales):
#         x = x + scale

#         # linear part
#         a1 = (et_2-et_1)/(em_2-em_1)
#         b1 = et_1-a1*em_1
#         linear = a1*x+b1

#         # exponential part
#         a2 = a1/a3
#         expo = a2*(np.exp(a3*(x-em_1))-1)+et_1

#         # adjust slightly
#         expo0 = a2*(np.exp(a3*(-em_1))-1)+et_1
#         linear = linear-expo0
#         expo = expo-expo0

#         results.append(np.where(x<em_1, expo, linear))
#     return np.concatenate(results)

def normalize_expr_df_to_rpm(df):
    # normalize
    df = df.div(df.sum(axis=0), axis=1) * 1000000
    # deduct the minimum and add one to the expression data to avoid division by 0
    df = df - df.min() + 1
    # normalize to rpm
    df = df.div(df.sum(axis=0), axis=1) * 1000000
    return df





