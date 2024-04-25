import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from .utils import ann_score_df
import scanpy as sc

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

# plt.rcParams['font.sans-serif'] = [
#     'Helvetica', 'Arial', 'Verdana', 'Bitstream Vera Sans'
# ]
# plt.rcParams['font.size'] = 8
# plt.rcParams['font.weight'] = 'regular'
# plt.rcParams['text.color'] = almost_black
#
# axisLineWidth = .5
# plt.rcParams['axes.linewidth'] = axisLineWidth
# plt.rcParams['lines.linewidth'] = 1.5
#
# plt.rcParams['axes.facecolor'] = 'white'
# plt.rcParams['axes.edgecolor'] = almost_black
# plt.rcParams['axes.labelcolor'] = almost_black
# # plt.rcParams['axes.color_cycle'] = dark2_all
#
# plt.rcParams['patch.edgecolor'] = 'none'
# plt.rcParams['patch.linewidth'] = .25
# # plt.rcParams['patch.facecolor'] = dark2_all[0]
#
# plt.rcParams['savefig.dpi'] = 1000
# plt.rcParams['savefig.format'] = 'svg'
#
# plt.rcParams['legend.frameon'] = False
# plt.rcParams['legend.handletextpad'] = .25
# plt.rcParams['legend.fontsize'] = 8
# plt.rcParams['legend.numpoints'] = 1
# plt.rcParams['legend.scatterpoints'] = 1
#
# plt.rcParams['ytick.direction'] = 'out'
# plt.rcParams['ytick.color'] = almost_black
# plt.rcParams['ytick.major.width'] = axisLineWidth
# plt.rcParams['xtick.direction'] = 'out'
# plt.rcParams['xtick.color'] = almost_black
# plt.rcParams['xtick.major.width'] = axisLineWidth


def draw_threshold(x, threshold, pseudo_sd):
    return threshold * pseudo_sd * (1 if x > 0 else -1) / abs(x)


def prep_data(df_in, threshold, ctrl_label):
    df = df_in.copy()

    df = ann_score_df(df, threshold=threshold, ctrl_label = ctrl_label)

    df['-log10(pvalue)'] = np.log10(df.pvalue) * -1

    return df


def plot_volcano(ax, df_in, threshold, up_hit='resistance_hit', down_hit='sensitivity_hit', 
                 ctrl_label = 'no-targeting',
                 dot_size=1,
                 xlim_l=-5, xlim_r=5,
                 ylim=6):
    df = prep_data(df_in, threshold, ctrl_label)

    # Scatter plot for each category
    ax.scatter(df.loc[df['label'] == 'target_non_hit', 'score'],
               df.loc[df['label'] == 'target_non_hit', '-log10(pvalue)'],
               alpha=0.1, s=dot_size, c='black', label='target_non_hit')
    ax.scatter(df.loc[df['label'] == up_hit, 'score'], df.loc[df['label'] == up_hit, '-log10(pvalue)'],
               alpha=0.9, s=dot_size, c='#fcae91', label=up_hit)
    ax.scatter(df.loc[df['label'] == down_hit, 'score'], df.loc[df['label'] == down_hit, '-log10(pvalue)'],
               alpha=0.9, s=dot_size, c='#bdd7e7', label=down_hit)
    ax.scatter(df.loc[df['label'] == 'non-targeting', 'score'],
               df.loc[df['label'] == 'non-targeting', '-log10(pvalue)'],
               alpha=0.1, s=dot_size, c='gray', label='non-targeting')

    # Set x-axis and y-axis labels
    ax.set_xlabel('phenotype score')
    ax.set_ylabel('-log10(p-value)')

    # Set x-axis limits
    ax.set_xlim(xlim_l, xlim_r)

    # Set y-axis limits
    ax.set_ylim(0.1, ylim)

    # Add legend
    ax.legend()


def label_as_black(ax, df_in, label, threshold, size=2, size_txt=None, 
                   ctrl_label = 'no-targeting',
                   t_x=.5, t_y=-0.1):
    df = prep_data(df_in, threshold, ctrl_label)

    target_data = df[df['target'] == label]

    # Scatter plot for labeled data
    ax.scatter(target_data['score'], target_data['-log10(pvalue)'],
               s=size, linewidth=0.5, edgecolors='black', facecolors='black', label='target')

    if not size_txt:
        size_txt = size * 2

    # Annotate the points
    for i, _ in enumerate(target_data['target']):
        txt = target_data['target'].iloc[i]
        ax.annotate(txt, (target_data['score'].iloc[i] + t_x, target_data['-log10(pvalue)'].iloc[i] + t_y),
                    color='black', size=size_txt)


def label_sensitivity_hit(ax, df_in, label, threshold, size=2, size_txt=None,
                          ctrl_label = 'no-targeting',
                          t_x=.5, t_y=-0.1):
    df = prep_data(df_in, threshold, ctrl_label)

    target_data = df[df['target'] == label]

    # Scatter plot for labeled data
    ax.scatter(target_data['score'], target_data['-log10(pvalue)'],
               s=size, linewidth=0.5, edgecolors='black', facecolors='#3182bd', label='target')

    if not size_txt:
        size_txt = size * 2

    # Annotate the points
    for i, _ in enumerate(target_data['target']):
        txt = target_data['target'].iloc[i]
        ax.annotate(txt, (target_data['score'].iloc[i] + t_x, target_data['-log10(pvalue)'].iloc[i] + t_y),
                    color='black', size=size_txt)


def label_resistance_hit(ax, df_in, label, threshold, size=2, size_txt=None, 
                         ctrl_label = 'no-targeting',
                         t_x=.5, t_y=-0.1):
    df = prep_data(df_in, threshold, ctrl_label)

    target_data = df[df['target'] == label]

    # Scatter plot for labeled data
    ax.scatter(target_data['score'], target_data['-log10(pvalue)'],
               s=size, linewidth=0.5, edgecolors='black', facecolors='#de2d26', label='target')

    if not size_txt:
        size_txt = size * 2

    # Annotate the points
    for i, _ in enumerate(target_data['target']):
        txt = target_data['target'].iloc[i]
        ax.annotate(txt, (target_data['score'].iloc[i] + t_x, target_data['-log10(pvalue)'].iloc[i] + t_y),
                    color='black', size=size_txt)


def plotReplicateScatter(ax, adat_in, x, y, title, min_val=None, max_val=None, log_transform=True):
    adat = adat_in[[x, y], :].copy()

    adat.obs.index = [f'Replicate {str(r)}' for r in adat.obs.replicate.to_list()]
    x_lab, y_lab = [f'Replicate {str(r)}' for r in adat.obs.replicate.to_list()]

    if log_transform:
        sc.pp.log1p(adat)
    
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
        ax=ax
    )
    ax.set_ylim(min_val, max_val)
    ax.set_xlim(min_val, max_val)
    ax.tick_params(axis='both', labelsize=10)
    ax.get_legend().remove()

    ax.grid(False)
