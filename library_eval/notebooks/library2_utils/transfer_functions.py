def transfer_function(x, c1, c2):
    """The expression is assumed to be normalized to one. The microRNA data is assumed to be linear.
    The return value is linear (not log10)"""
    c1 = 10**c1
    c2 = 10**c2
    result = (1 / (1 + x / c1)) * (1 + x / c2)
    if type(result) == float:
        return result
    else:
        return result.astype(float)

def inverse_transfer(expr, c1, c2):
    """Inverts the transfer function. Expects expression values between 0 and 1."""
    c1 = 10**c1
    c2 = 10**c2
    
    result = c1*c2*(expr - 1) / (c1 - c2*expr)
    return result.astype(float)