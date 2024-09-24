import numpy as np

def tsi(x):
    """This function calculates the TSI for a given combination of cell lines.
    The input is a numpy array with the expression of the cell lines in the columns
    and the microRNAs in the rows."""
    # check if the orientation of x is correct
    if x.shape[1] > 50:
        raise ValueError("The number of columns is very large. The columns are supposed to be cell lines. Is this the case?")

    # if x is not normalized yet, normalize it
    # THIS IS IMPORTANT
    x = x/x.max(axis=1, keepdims=True)

    tsi = np.sum(1-x, axis=1)/(x.shape[1]-1)
    
    return tsi

def calculate_quality(df, cell_line):
    """This functions calculates the quality of a design for a given cell line."""
    quality = df[cell_line]*df["tsi"]

    return quality