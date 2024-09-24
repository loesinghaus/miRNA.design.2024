def transfer_function(x, c1=10**3.6, c2=10**10):
    """The expression is assumed to be normalized to one. The microRNA data is assumed to be linear."""
    result = (1 / (1 + x / c1)) * (1 + x / c2)
    return result

def inverse_transfer(expr, c1=10**3.6, c2=10**10):
    return c1*c2*(expr - 1) / (c1 - c2*expr)