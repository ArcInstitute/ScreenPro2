import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ._utils import yellow_blue
from ._rank import rank_plot


def volcano_plot(
            ax, df, up_hit, down_hit,
            score_col='score', pvalue_col='pvalue',
            xlabel='phenotype score',
            ylabel='-log10(pvalue)',
            dot_size=1,
            xlims='auto',
            ylims='auto',
            ctrl_label='negative_control',
            **args
            ):
        
        if f'-log10({pvalue_col})' not in df.columns:
            df[f'-log10({pvalue_col})'] = np.log10(df[pvalue_col]) * -1

        if xlims == 'auto':
            xlims = (df[score_col].min() - 0.1, df[score_col].max() + 0.1)
        if ylims == 'auto':
            ylims = (df[f'-log10({pvalue_col})'].min() - 0.1, df[f'-log10({pvalue_col})'].max() + 0.1)
        
        # Scatter plot for each category
        ax.scatter( df.loc[df['label'] == 'target_non_hit', score_col],
                    df.loc[df['label'] == 'target_non_hit', f'-log10({pvalue_col})'],
                    alpha=0.1, s=dot_size, c='black', label='target_non_hit',
                    **args)
        
        ax.scatter( df.loc[df['label'] == up_hit, score_col], 
                    df.loc[df['label'] == up_hit, f'-log10({pvalue_col})'],
                    alpha=0.9, s=dot_size, c='#fcae91', label=up_hit,
                    **args)
        
        ax.scatter( df.loc[df['label'] == down_hit, score_col], 
                    df.loc[df['label'] == down_hit, f'-log10({pvalue_col})'],
                    alpha=0.9, s=dot_size, c='#bdd7e7', label=down_hit,
                    **args)
        
        ax.scatter( df.loc[df['label'] == ctrl_label, score_col],
                    df.loc[df['label'] == ctrl_label, f'-log10({pvalue_col})'],
                    alpha=0.1, s=dot_size, c='gray', label=ctrl_label,
                    **args)

        # Set x-axis and y-axis labels
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        # Set x-axis limits
        ax.set_xlim(xlims)

        # Set y-axis limits
        ax.set_ylim(ylims)

        # Add legend
        ax.legend()


def label_by_color(ax, df_in, label,
                   x_col, y_col,
                   size=2, size_txt="auto",
                   edgecolors='black', facecolors='black',
                   textcolor='black',
                   t_x=.5, t_y=-0.1, **args):

    if size_txt == 'auto': size_txt = size * 2

    df = df_in.copy()
    target_data = df[df['target'] == label]

    # Scatter plot for labeled data
    ax.scatter(
        target_data[x_col], target_data[y_col],
        s=size, linewidth=0.5, 
        edgecolors=edgecolors, 
        facecolors=facecolors, label='target',
        **args
    )

    if size_txt == None:
        pass
    else:
        # Annotate the points
        for i, _ in enumerate(target_data['target']):
            txt = target_data['target'].iloc[i]
            ax.annotate(txt, (target_data[x_col].iloc[i] + t_x, target_data[y_col].iloc[i] + t_y),
                        color=textcolor, size=size_txt)


def label_as_black(ax, df_in, label, 
                    x_col='score', y_col='-log10(pvalue)',
                    size=2, size_txt="auto",
                    t_x=.5, t_y=-0.1,
                    **args):
    label_by_color(
        ax, df_in, label, 
        x_col=x_col, y_col=y_col,
        size=size, size_txt=size_txt,
        edgecolors='black', facecolors='black',
        textcolor='black',
        t_x=t_x, t_y=t_y,
        **args)


def label_sensitivity_hit(ax, df_in, label, 
                            x_col='score', y_col='-log10(pvalue)',
                            size=2, size_txt="auto",
                            t_x=.5, t_y=-0.1,
                            **args):
    label_by_color(
        ax, df_in, label, 
        x_col=x_col, y_col=y_col,
        size=size, size_txt=size_txt,
        edgecolors='black', facecolors='#3182bd',
        textcolor='black',
        t_x=t_x, t_y=t_y,
        **args)


def label_resistance_hit(ax, df_in, label, 
                            x_col='score', y_col='-log10(pvalue)',
                            size=2, size_txt="auto",
                            t_x=.5, t_y=-0.1,
                            **args):
    label_by_color(
        ax, df_in, label, 
        x_col=x_col, y_col=y_col,
        size=size, size_txt=size_txt,
        edgecolors='black', facecolors='#de2d26',
        textcolor='black',
        t_x=t_x, t_y=t_y,
        **args)


# def drawRhoGammaScatter(
#         self, ax,
#         rho_df=None, gamma_df=None,
#         dot_size=1,
#         score_col='score',
#         xlabel='auto',
#         ylabel='auto',
#         xlims='auto',
#         ylims='auto',
#         **args
#         ):
#     #TODO: fix by making a single dataframe with both rho and gamma scores
#     if rho_df is None:
#         _, _, rho_df = self._prep_data()
#     if gamma_df is None:
#         gamma_df, _, _ = self._prep_data()

#     if xlabel == 'auto':
#         xlabel = self.rho_score_name.replace(':', ': ').replace('_', ' ')
#     if ylabel == 'auto':
#         ylabel = self.gamma_score_name.replace(':', ': ').replace('_', ' ')
    
#     # color by rho score labels
#     up_hit = 'resistance_hit'
#     down_hit = 'sensitivity_hit'

#     # Scatter plot for each category
#     ax.scatter( rho_df.loc[rho_df['label'] == 'target_non_hit', score_col],
#                 gamma_df.loc[rho_df['label'] == 'target_non_hit', score_col],
#                 alpha=0.1, s=dot_size, c='black', label='target_non_hit',
#                 **args)
    
#     ax.scatter( rho_df.loc[rho_df['label'] == up_hit, score_col],
#                 gamma_df.loc[rho_df['label'] == up_hit, score_col],
#                 alpha=0.9, s=dot_size, c='#fcae91', label=up_hit,
#                 **args)
    
#     ax.scatter( rho_df.loc[rho_df['label'] == down_hit, score_col],
#                 gamma_df.loc[rho_df['label'] == down_hit, score_col],
#                 alpha=0.9, s=dot_size, c='#bdd7e7', label=down_hit,
#                 **args)
    
#     ax.scatter( rho_df.loc[rho_df['label'] == self.ctrl_label, score_col],
#                 gamma_df.loc[rho_df['label'] == self.ctrl_label, score_col],
#                 alpha=0.1, s=dot_size, c='gray', label=self.ctrl_label,
#                 **args)
    
#     # Set x-axis and y-axis labels
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel(ylabel)

#     # Set x-axis limits
#     ax.set_xlim(xlims)
#     ax.set_ylim(ylims)

#     # Add legend
#     ax.legend()
    
