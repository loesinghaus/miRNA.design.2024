import pandas as pd
import numpy as np

def normalize_expr_df_to_rpm(df, minimum="none"):
    # is it log
    log = False
    if df.max().max() < 10:
        df = 10**df
        log = True

    if minimum == "none":
        minimum = 1
    elif log:
        minimum = 10**minimum

    # normalize
    df = df.div(df.sum(axis=0), axis=1) * 1000000
    # deduct the minimum and add one to the expression data to avoid division by 0
    df = df - df.min() + minimum
    # normalize to rpm
    df = df.div(df.sum(axis=0), axis=1) * 1000000

    if log:
        df = np.log10(df)

    return df


def normalize_expr_df_to_rpm_with_index(df, index):
    # is it log
    log = False
    if df.max().max() < 10:
        df = 10**df
        log = True
    
    # normalize
    df_subset = df[df.index.isin(index)].copy()
    
    df = df.div(df_subset.sum(axis=0), axis=1) * 1000000
    df_subset = df_subset.div(df_subset.sum(axis=0), axis=1) * 1000000
    
    # deduct the minimum and add one to the expression data to avoid division by 0
    df = df - df_subset.min().min() + 1
    df_subset = df_subset - df_subset.min().min() + 1
    
    # normalize to rpm
    df = df.div(df_subset.sum(axis=0), axis=1) * 1000000
    df_subset = df_subset.div(df_subset.sum(axis=0), axis=1) * 1000000
    
    if log:
        df = np.log10(df)
    
    return df


def normalize_expr_df_to_rpm_with_partner(df_partner, df, minimum="none"):
    # is it log
    log = False
    if df.max().max() < 10:
        df = 10**df
        df_partner = 10**df_partner
        log = True
    
    if minimum == "none":
        minimum = 1
    elif log:
        minimum = 10**minimum
    
    df = df.div(df_partner.sum(axis=0), axis=1) * 1000000
    df_partner = df_partner.div(df_partner.sum(axis=0), axis=1) * 1000000
    
    # deduct the minimum and add one to the expression data to avoid division by 0
    df = df - df_partner.min().min() + minimum
    df_partner = df_partner - df_partner.min().min() + minimum
    
    # normalize to rpm
    df = df.div(df_partner.sum(axis=0), axis=1) * 1000000
    df_partner = df_partner.div(df_partner.sum(axis=0), axis=1) * 1000000
    
    if log:
        df = np.log10(df)
        df_partner = np.log10(df_partner)
    
    return df_partner, df


