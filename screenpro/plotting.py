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


## Scatter plot of replicates

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


## Phenotype plotting functions

class PheScatterPlot:

    def __init__(self,threshold=3,ctrl_label='no-targeting') -> None:
        self.threshold = threshold
        self.ctrl_label = ctrl_label

    def _prepare_score_df(self, df_in):
        df = df_in.copy()

        df = ann_score_df(df, self.threshold, self.ctrl_label)

        df['-log10(pvalue)'] = np.log10(df.pvalue) * -1

        return df
    
    def _label_by_color(self, ax, df_in, label, size=2, size_txt="auto",
                        x_col='score', y_col='-log10(pvalue)',
                        edgecolors='black', facecolors='black',
                        textcolor='black',
                        t_x=.5, t_y=-0.1):
        df = self._prepare_score_df(df_in, self.threshold, self.ctrl_label)

        target_data = df[df['target'] == label]

        # Scatter plot for labeled data
        ax.scatter(
            target_data[x_col], target_data[y_col],
            s=size, linewidth=0.5, 
            edgecolors=edgecolors, 
            facecolors=facecolors, label='target'
        )

        if size_txt == None:
            pass
        elif size_txt == 'auto':
            size_txt = size * 2
        else:
            # Annotate the points
            for i, _ in enumerate(target_data['target']):
                txt = target_data['target'].iloc[i]
                ax.annotate(txt, (target_data[x_col].iloc[i] + t_x, target_data[y_col].iloc[i] + t_y),
                            color=textcolor, size=size_txt)


class Volcano(PheScatterPlot):

    def __init__(self, threshold=3, ctrl_label='no-targeting') -> None:
        super().__init__(threshold, ctrl_label)

    def plot(self, ax, df_in, up_hit='resistance_hit', down_hit='sensitivity_hit', 
             xlabel='phenotype score',
             ylabel='-log10(pvalue)',
             dot_size=1,
             xlims=(-5, 5),
             ylim=6):
        df = self._prepare_score_df(df_in)
        xlim_l, xlim_r = xlims

        # Scatter plot for each category
        ax.scatter( df.loc[df['label'] == 'target_non_hit', 'score'],
                    df.loc[df['label'] == 'target_non_hit', '-log10(pvalue)'],
                    alpha=0.1, s=dot_size, c='black', label='target_non_hit')
        
        ax.scatter( df.loc[df['label'] == up_hit, 'score'], 
                    df.loc[df['label'] == up_hit, '-log10(pvalue)'],
                    alpha=0.9, s=dot_size, c='#fcae91', label=up_hit)
        
        ax.scatter( df.loc[df['label'] == down_hit, 'score'], 
                    df.loc[df['label'] == down_hit, '-log10(pvalue)'],
                    alpha=0.9, s=dot_size, c='#bdd7e7', label=down_hit)
        
        ax.scatter( df.loc[df['label'] == self.ctrl_label, 'score'],
                    df.loc[df['label'] == self.ctrl_label, '-log10(pvalue)'],
                    alpha=0.1, s=dot_size, c='gray', label=self.ctrl_label)

        # Set x-axis and y-axis labels
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        # Set x-axis limits
        ax.set_xlim(xlim_l, xlim_r)

        # Set y-axis limits
        ax.set_ylim(0.1, ylim)

        # Add legend
        ax.legend()

    def label_as_black(self, ax, df_in, label, size=2, size_txt="auto",
                       t_x=.5, t_y=-0.1):
        self._label_by_color(ax, df_in, label, size, size_txt,
                             edgecolors='black', facecolors='black',
                             textcolor='black',
                             t_x=t_x, t_y=t_y)
    
    def label_sensitivity_hit(self, ax, df_in, label, size=2, size_txt="auto",
                              t_x=.5, t_y=-0.1):
        self._label_by_color(ax, df_in, label, size, size_txt,
                             edgecolors='black', facecolors='#3182bd',
                             textcolor='black',
                             t_x=t_x, t_y=t_y)
    
    def label_resistance_hit(self, ax, df_in, label, size=2, size_txt="auto",
                             t_x=.5, t_y=-0.1):
        self._label_by_color(ax, df_in, label, size, size_txt,
                             edgecolors='black', facecolors='#de2d26',
                             textcolor='black',
                             t_x=t_x, t_y=t_y)


class RhoGammaScatter(PheScatterPlot):
    def __init__(self, screen, threshold=3, ctrl_label='no-targeting') -> None:
        self.screen = screen
        super().__init__(threshold, ctrl_label)

    def _prep_data(self):

        gamma = self.screen.getPhenotypeScores(
            run_name=self.run_name,
            score_name=self.gamma_score_name,
            ctrl_label=self.ctrl_label,
            threshold=self.threshold,
        )

        rho = self.screen.getPhenotypeScores(
            run_name=self.run_name,
            score_name=self.rho_score_name,
            ctrl_label=self.ctrl_label,
            threshold=self.threshold
        )

        rho = self._prepare_score_df(rho)
        gamma = self._prepare_score_df(gamma)

        return rho, gamma
        
    def plot(self, ax,
             up_hit='resistance_hit', down_hit='sensitivity_hit',
             dot_size=1,
             xlabel=None,
             ylabel=None,
             xlims=(-5, 5),
             ylims=(-5, 5)):

        rho_df, gamma_df = self._prep_data(self.screen)

        # Scatter plot for each category
        ax.scatter( rho_df.loc[rho_df['label'] == 'target_non_hit', 'score'],
                    gamma_df.loc[rho_df['label'] == 'target_non_hit', 'score'],
                    alpha=0.1, s=dot_size, c='black', label='target_non_hit')
        
        ax.scatter( rho_df.loc[rho_df['label'] == up_hit, 'score'],
                    gamma_df.loc[rho_df['label'] == up_hit, 'score'],
                    alpha=0.9, s=dot_size, c='#fcae91', label=up_hit)
        
        ax.scatter( rho_df.loc[rho_df['label'] == down_hit, 'score'],
                    gamma_df.loc[rho_df['label'] == down_hit, 'score'],
                    alpha=0.9, s=dot_size, c='#bdd7e7', label=down_hit)
        
        ax.scatter( rho_df.loc[rho_df['label'] == self.ctrl_label, 'score'],
                    gamma_df.loc[rho_df['label'] == self.ctrl_label, 'score'],
                    alpha=0.1, s=dot_size, c='gray', label=self.ctrl_label)
        
        # Set x-axis and y-axis labels
        if not xlabel:
            ax.set_xlabel('rho score')
        else:
            ax.set_xlabel(xlabel)
        if not ylabel:
            ax.set_ylabel('gamma score')
        else:
            ax.set_ylabel(ylabel)

        # Set x-axis limits
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)

        # Add legend
        ax.legend()
    
    def label_as_black(self, ax, label, size=2, size_txt="auto",
                       t_x=.5, t_y=-0.1):
        self._label_by_color(ax, label, size, size_txt,
                             edgecolors='black', facecolors='black',
                             textcolor='black',
                             t_x=t_x, t_y=t_y)
        
    def label_sensitivity_hit(self, ax, label, size=2, size_txt="auto",
                              t_x=.5, t_y=-0.1):
    
        self._label_by_color(ax, label, size, size_txt,
                             edgecolors='black', facecolors='#3182bd',
                             textcolor='black',
                             t_x=t_x, t_y=t_y)
        
    def label_resistance_hit(self, ax, label, size=2, size_txt="auto",
                             t_x=.5, t_y=-0.1):
        
        self._label_by_color(ax, label, size, size_txt,
                             edgecolors='black', facecolors='#de2d26',
                             textcolor='black',
                             t_x=t_x, t_y=t_y)
