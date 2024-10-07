from matplotlib.legend_handler import HandlerPathCollection
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import numpy as np

class HandlerSize(HandlerPathCollection):
    def __init__(self, marker_size=12):
        self.marker_size = marker_size
        super().__init__()
        
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        markers = super().create_artists(legend, orig_handle,
                                         xdescent, ydescent, width, height, fontsize, trans)
        for marker in markers:
            marker.set_sizes([self.marker_size])  # Set desired legend marker size here
        return markers

def density_scatter(x, y, ax=None, sort=True, bins=20, **kwargs):
    """
    Scatter plot colored by 2d histogram
    """
    if ax is None:
        fig, ax = plt.subplots()
    
    data, x_e, y_e = np.histogram2d(x, y, bins=bins)
    z = gaussian_kde(np.vstack([x, y]))(np.vstack([x, y]))
    
    # Sort the points by density, so that the densest points are plotted last
    if sort:
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
    
    sc = ax.scatter(x, y, c=z, **kwargs)
    
    return ax, sc
    
def return_pvalue_text(p_value):
    formatted_p_value = f"{p_value:.2e}"
    base, exponent = formatted_p_value.split('e')
    exponent = int(exponent)
    return f'$p = {base} \\cdot 10^{{\\mathrm{{{exponent}}}}}$'