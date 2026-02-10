import numpy as np

def transfer_function_linear(x, c1):
    """The expression is assumed to be normalized to one. The microRNA data is assumed to be linear.
    The return value is linear (not log10). c1 is assumed to be in log10."""
    c1 = 10**c1
    result = (1 / (1 + x / c1)) 
    if type(result) == float:
        return result
    else:
        return result.astype(float)
    
def transfer_function_log(x, c1):
    """The expression is assumed to be normalized to one. The microRNA data is assumed to be log10.
    The return value is log10. c1 is assumed to be in log10."""
    c1 = 10**c1
    x = 10**x
    result = np.log10(1 / (1 + x / c1)) 
    if type(result) == float:
        return result
    else:
        return result.astype(float)