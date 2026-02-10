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