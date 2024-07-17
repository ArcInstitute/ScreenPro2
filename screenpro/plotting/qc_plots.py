import numpy as np
import scanpy as sc
from ._utils import almost_black


## Histogram of guide counts distribution
def plotCountDistribution(ax, adat, title, **args):
    pass


## Scatter plot of replicates
def plotReplicateScatter(ax, adat_in, x, y, title, min_val=None, max_val=None, log_transform=True, **args):
    adat = adat_in[[x, y], :].copy()

    adat.obs.index = [f'Replicate {str(r)}' for r in adat.obs.replicate.to_list()]
    x_lab, y_lab = [f'Replicate {str(r)}' for r in adat.obs.replicate.to_list()]

    if log_transform:
        adat.X = np.log10(adat.X+1)
    
    if min_val is None:
        min_val = min([adat.to_df().loc[x_lab,:].min(), adat.to_df().loc[y_lab,:].min()])
        min_val = min_val * 1.1 
    if max_val is None:
        max_val = max([adat.to_df().loc[x_lab,:].max(), adat.to_df().loc[y_lab,:].max()])
        max_val = max_val * 1.1
    
    sc.pl.scatter(
        adat,
        x_lab, y_lab,
        legend_fontsize='xx-large',
        palette=[almost_black, '#BFBFBF'],
        color='targetType',
        title=title,
        size=5,
        show=False,
        ax=ax,
        **args
    )
    ax.set_ylim(min_val, max_val)
    ax.set_xlim(min_val, max_val)
    ax.tick_params(axis='both', labelsize=10)
    ax.get_legend().remove()

    ax.grid(False)
