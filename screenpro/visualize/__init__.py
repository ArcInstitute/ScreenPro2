## Copyright (c) 2022-2024 ScreenPro2 Development Team.
## All rights reserved.
## Gilbart Lab, UCSF / Arc Institute.
## Multi-Omics Tech Center, Arc Insititue.

import numpy as np
import scanpy as sc
from .qc_plots import *

## Phenotype plotting functions
class DrugScreenPlotter:
    def __init__(self, screen, treated, untreated, t0='T0', threshold=3, ctrl_label='negative_control',run_name='auto'):
        self.screen = screen
        self.threshold = threshold
        self.ctrl_label = ctrl_label
        self.run_name = run_name
        self.gamma_score_name = f'gamma:{untreated}_vs_{t0}'
        self.rho_score_name = f'rho:{treated}_vs_{untreated}'
        self.tau_score_name = f'tau:{treated}_vs_{t0}'

    def _label_by_color(self, ax, df_in, label, 
                        x_col, y_col,
                        size=2, size_txt="auto",
                        edgecolors='black', facecolors='black',
                        textcolor='black',
                        t_x=.5, t_y=-0.1, **args):
        
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
        elif size_txt == 'auto':
            size_txt = size * 2
        else:
            # Annotate the points
            for i, _ in enumerate(target_data['target']):
                txt = target_data['target'].iloc[i]
                ax.annotate(txt, (target_data[x_col].iloc[i] + t_x, target_data[y_col].iloc[i] + t_y),
                            color=textcolor, size=size_txt)

    def _prep_data(self, score_col='score', pvalue_col='pvalue'):

        gamma = self.screen.getPhenotypeScores(
            run_name=self.run_name,
            score_name=self.gamma_score_name,
            threshold=self.threshold,
            ctrl_label=self.ctrl_label,
            score_col=score_col,
            pvalue_col=pvalue_col
        )

        gamma[f'-log10({pvalue_col})'] = np.log10(gamma[pvalue_col]) * -1

        tau = self.screen.getPhenotypeScores(
            run_name=self.run_name,
            score_name=self.tau_score_name,
            threshold=self.threshold,
            ctrl_label=self.ctrl_label,
            score_col=score_col,
            pvalue_col=pvalue_col
        )
        tau[f'-log10({pvalue_col})'] = np.log10(tau[pvalue_col]) * -1

        rho = self.screen.getPhenotypeScores(
            run_name=self.run_name,
            score_name=self.rho_score_name,
            threshold=self.threshold,
            ctrl_label=self.ctrl_label,
            score_col=score_col,
            pvalue_col=pvalue_col
        )
        rho[f'-log10({pvalue_col})'] = np.log10(rho[pvalue_col]) * -1

        return gamma, tau, rho
    
    def _volcano(
            self, ax, df, up_hit, down_hit,
            xlabel='phenotype score',
            ylabel='-log10(pvalue)',
            dot_size=1,
            xlims=(-5, 5),
            ylims=(0.1,6),
            **args
            ):
        
        # Scatter plot for each category
        ax.scatter( df.loc[df['label'] == 'target_non_hit', 'score'],
                    df.loc[df['label'] == 'target_non_hit', '-log10(pvalue)'],
                    alpha=0.1, s=dot_size, c='black', label='target_non_hit',
                    **args)
        
        ax.scatter( df.loc[df['label'] == up_hit, 'score'], 
                    df.loc[df['label'] == up_hit, '-log10(pvalue)'],
                    alpha=0.9, s=dot_size, c='#fcae91', label=up_hit,
                    **args)
        
        ax.scatter( df.loc[df['label'] == down_hit, 'score'], 
                    df.loc[df['label'] == down_hit, '-log10(pvalue)'],
                    alpha=0.9, s=dot_size, c='#bdd7e7', label=down_hit,
                    **args)
        
        ax.scatter( df.loc[df['label'] == self.ctrl_label, 'score'],
                    df.loc[df['label'] == self.ctrl_label, '-log10(pvalue)'],
                    alpha=0.1, s=dot_size, c='gray', label=self.ctrl_label,
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

    def drawVolcanoRho(
            self, ax,
            rho_df=None,
            dot_size=1,
            xlabel='auto',
            ylabel='-log10(pvalue)',
            xlims=(-5, 5),
            ylims=(0.1, 6),
            **args
            ):
        if rho_df is None:
            _, _, rho_df = self._prep_data()
        if xlabel == 'auto':
            xlabel = self.rho_score_name.replace(':', ': ').replace('_', ' ')
        
        self._volcano(ax, rho_df, 
                      up_hit='resistance_hit', down_hit='sensitivity_hit',
                      xlabel=xlabel, ylabel=ylabel,
                      dot_size=dot_size, xlims=xlims, ylims=ylims,
                      **args)
    
    def drawVolcanoGamma(
            self, ax,
            gamma_df=None,
            dot_size=1,
            xlabel='auto',
            ylabel='-log10(pvalue)',
            xlims=(-5, 5),
            ylims=(0.1, 6),
            **args
            ):
        if gamma_df is None:
            gamma_df, _, _ = self._prep_data()
        if xlabel == 'auto':
            xlabel = self.gamma_score_name.replace(':', ': ').replace('_', ' ')
        
        self._volcano(ax, gamma_df, 
                      up_hit='up_hit', down_hit='essential_hit',
                      xlabel=xlabel, ylabel=ylabel,
                      dot_size=dot_size, xlims=xlims, ylims=ylims,
                      **args)
        
    def drawVolcanoTau(
            self, ax,
            tau_df=None,
            dot_size=1,
            xlabel='auto',
            ylabel='-log10(pvalue)',
            xlims=(-5, 5),
            ylims=(0.1, 6),
            **args
            ):
        if tau_df is None:
            _, tau_df, _, = self._prep_data()
        if xlabel == 'auto':
            xlabel = self.tau_score_name.replace(':', ': ').replace('_', ' ')
        
        self._volcano(ax, tau_df, 
                      up_hit='up_hit', down_hit='down_hit',
                      xlabel=xlabel, ylabel=ylabel,
                      dot_size=dot_size, xlims=xlims, ylims=ylims,
                      **args)

    def drawRhoGammaScatter(
            self, ax,
            rho_df=None, gamma_df=None,
            dot_size=1,
            xlabel='auto',
            ylabel='auto',
            xlims=(-5, 5),
            ylims=(-5, 5),
            **args
            ):
        #TODO: fix by making a single dataframe with both rho and gamma scores
        if rho_df is None:
            _, _, rho_df = self._prep_data()
        if gamma_df is None:
            gamma_df, _, _ = self._prep_data()

        if xlabel == 'auto':
            xlabel = self.rho_score_name.replace(':', ': ').replace('_', ' ')
        if ylabel == 'auto':
            ylabel = self.gamma_score_name.replace(':', ': ').replace('_', ' ')
        
        # color by rho score labels
        up_hit = 'resistance_hit'
        down_hit = 'sensitivity_hit'

        # Scatter plot for each category
        ax.scatter( rho_df.loc[rho_df['label'] == 'target_non_hit', 'score'],
                    gamma_df.loc[rho_df['label'] == 'target_non_hit', 'score'],
                    alpha=0.1, s=dot_size, c='black', label='target_non_hit',
                    **args)
        
        ax.scatter( rho_df.loc[rho_df['label'] == up_hit, 'score'],
                    gamma_df.loc[rho_df['label'] == up_hit, 'score'],
                    alpha=0.9, s=dot_size, c='#fcae91', label=up_hit,
                    **args)
        
        ax.scatter( rho_df.loc[rho_df['label'] == down_hit, 'score'],
                    gamma_df.loc[rho_df['label'] == down_hit, 'score'],
                    alpha=0.9, s=dot_size, c='#bdd7e7', label=down_hit,
                    **args)
        
        ax.scatter( rho_df.loc[rho_df['label'] == self.ctrl_label, 'score'],
                    gamma_df.loc[rho_df['label'] == self.ctrl_label, 'score'],
                    alpha=0.1, s=dot_size, c='gray', label=self.ctrl_label,
                    **args)
        
        # Set x-axis and y-axis labels
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        # Set x-axis limits
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)

        # Add legend
        ax.legend()
        
    def label_as_black(self, ax, df_in, label, 
                       x_col='score', y_col='-log10(pvalue)',
                       size=2, size_txt="auto",
                       t_x=.5, t_y=-0.1,
                       **args):
        self._label_by_color(ax, df_in, label, 
                             x_col=x_col, y_col=y_col,
                             size=size, size_txt=size_txt,
                             edgecolors='black', facecolors='black',
                             textcolor='black',
                             t_x=t_x, t_y=t_y,
                             **args)
    
    def label_sensitivity_hit(self, ax, df_in, label, 
                              x_col='score', y_col='-log10(pvalue)',
                              size=2, size_txt="auto",
                              t_x=.5, t_y=-0.1,
                              **args):
        self._label_by_color(ax, df_in, label, 
                             x_col=x_col, y_col=y_col,
                             size=size, size_txt=size_txt,
                             edgecolors='black', facecolors='#3182bd',
                             textcolor='black',
                             t_x=t_x, t_y=t_y,
                             **args)
    
    def label_resistance_hit(self, ax, df_in, label, 
                             x_col='score', y_col='-log10(pvalue)',
                             size=2, size_txt="auto",
                             t_x=.5, t_y=-0.1,
                             **args):
        self._label_by_color(ax, df_in, label, 
                             x_col=x_col, y_col=y_col,
                             size=size, size_txt=size_txt,
                             edgecolors='black', facecolors='#de2d26',
                             textcolor='black',
                             t_x=t_x, t_y=t_y,
                             **args)
