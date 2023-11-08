# import
import pandas as pd
import numpy as np
from plotnine import *          # <-- NOT recommended!
import matplotlib.pyplot as plt
import matplotlib

# variables
almost_black = '#111111'
dark2 = ['#1b9e77',
         '#d95f02',
         '#7570b3',
         '#e7298a',
         '#66a61e',
         '#e6ab02',
         '#a6761d',
         '#666666']
blue_yellow = matplotlib.colors.LinearSegmentedColormap.from_list(
    'BuYl', [(0, '#ffff00'), (.49, '#000000'), (.51, '#000000'), (1, '#0000ff')])
blue_yellow.set_bad('#999999', 1)
yellow_blue = matplotlib.colors.LinearSegmentedColormap.from_list(
    'YlBu', [(0, '#0000ff'), (.49, '#000000'), (.51, '#000000'), (1, '#ffff00')])
yellow_blue.set_bad('#999999', 1)

plt.rcParams['font.sans-serif'] = [
    'Helvetica', 'Arial', 'Verdana', 'Bitstream Vera Sans'
]
plt.rcParams['font.size'] = 8
plt.rcParams['font.weight'] = 'regular'
plt.rcParams['text.color'] = almost_black

axisLineWidth = .5
plt.rcParams['axes.linewidth'] = axisLineWidth
plt.rcParams['lines.linewidth'] = 1.5

plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['axes.edgecolor'] = almost_black
plt.rcParams['axes.labelcolor'] = almost_black
# plt.rcParams['axes.color_cycle'] = dark2_all

plt.rcParams['patch.edgecolor'] = 'none'
plt.rcParams['patch.linewidth'] = .25
# plt.rcParams['patch.facecolor'] = dark2_all[0]

plt.rcParams['savefig.dpi'] = 1000
plt.rcParams['savefig.format'] = 'svg'

plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.handletextpad'] = .25
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['legend.numpoints'] = 1
plt.rcParams['legend.scatterpoints'] = 1

plt.rcParams['ytick.direction'] = 'out'
plt.rcParams['ytick.color'] = almost_black
plt.rcParams['ytick.major.width'] = axisLineWidth
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['xtick.color'] = almost_black
plt.rcParams['xtick.major.width'] = axisLineWidth


def plot_ggplot_scatter(adata, x, y, color, size, alpha, shape, facet):
    """
    Plot scatter plot using ggplot2 style via plotnine.
    Args:
        adata: anndata object
        x: x-axis variable
        y: y-axis variable
        color: color variable
        size: size variable
        alpha: alpha variable
        shape: shape variable
        facet: facet variable
    """
    scatter_p = (
        ggplot(adata.obs)
            + geom_point(
                aes(
                    x = x,
                    y = y,
                    color = color,
                    size = size,
                    alpha = alpha,
                    shape = shape
                    ), 
                stroke=0.2
            )
            + facet_wrap(facet)
            + theme_classic()
            + theme(
                panel_grid_major   = element_blank(),
                 panel_grid_minor   = element_blank(),
                panel_background   = element_blank(),

                legend_background  = element_blank(),
                legend_position    = 'top',
                legend_direction   = 'horizontal', # affected by the ncol=2
                legend_text_legend = element_text(size=8),

                axis_line          = element_line(size=2),
                axis_text_x        = element_text(size=8),
                axis_text_y        = element_text(size=8),
                axis_title_x       = element_text(weight='bold', size=12),
                axis_title_y       = element_text(weight='bold', size=12),

                text               = element_text(font = 'arial'),

                figure_size        = (4.5, 5)
            )
            + xlim(-60,80)
    )

    return scatter_p


def plot_ggplot_pca(adata):
    """
    Plot PCA using ggplot2 style via plotnine.
    Args:
        adata: anndata object
    """
    # Create a dataframe with the PCA coordinates and the metadata
    pca = pd.concat([
        pd.DataFrame(
            adata.obsm['X_pca'][:,[0,1]],
            index=adata.obs.index,
            columns=['PC-1', 'PC-2']
        ),
        adata.obs.drop('replicate', axis=1)
    ], axis=1)

    pca_p = (
        ggplot(pca)
            + geom_point(
                aes(
                    x = 'PC-1',
                    y = 'PC-2',
                    fill='score',
                    # shape = 'treatment'
                ), 
                color='black', 
                size=8
            )
            + geom_text(
                aes(
                    x='PC-1',
                    y='PC-2',
                    label='score',
                    size=3
                ),
                nudge_y=3,
                nudge_x=7,
            )
            + theme_classic()
            + theme(
                panel_grid_major   = element_blank(),
                panel_grid_minor   = element_blank(),
                panel_background   = element_blank(),

                legend_background  = element_blank(),
                legend_position    = 'top',
                legend_direction   = 'horizontal', # affected by the ncol=2
                legend_text_legend = element_text(size=8),

                axis_line          = element_line(size=2),
                axis_text_x        = element_text(size=8),
                axis_text_y        = element_text(size=8),
                axis_title_x       = element_text(weight='bold', size=12),
                axis_title_y       = element_text(weight='bold', size=12),

                text               = element_text(font = 'arial'),

                figure_size        = (4.5, 5)
            )
            + xlim(-60,80)
    )

    return pca_p