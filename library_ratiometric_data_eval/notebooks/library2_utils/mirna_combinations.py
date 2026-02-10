def match_design_name(x, design_miRNA_dict):
    """Takes a list of miRNA names (x) and returns the design name in design_miRNA_dict.
    Designed to be used with the apply function on a dataframe."""
    
    # sort the miRNA names alphabetically
    x = sorted(x)
    # make them a tuple
    x = tuple(x)
    # if the tuple is in the design_miRNA_dict, return the design name
    if x in design_miRNA_dict.keys():
        return design_miRNA_dict[x]
    else:
        return None

def get_combinations(df):
    """This function takes a dataframe with miRNA name columns and returns a list of tuples."""
    
    # get all columns that contain "miRNA"
    miRNA_columns = [column for column in df.columns if "miRNA" in column]
    # make a list of tuples of the miRNA names
    miRNA_combinations = list(zip(*[df[miRNA_column] for miRNA_column in miRNA_columns]))
    
    return miRNA_combinations