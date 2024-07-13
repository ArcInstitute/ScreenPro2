import numpy as np
import pandas as pd
import bokeh
import bokeh.plotting


class DataDashboard:

    def __init__(self):
        pass

    def _new_plot(self,title,tooltips,width,height,toolbar_location):
        
        TOOLS = "box_select,box_zoom,lasso_select,reset,save,wheel_zoom,pan,copy,undo,redo,reset,examine,fullscreen"

        # create a new plot with a specific size
        p = bokeh.plotting.figure(
            sizing_mode="stretch_width",
            tools=TOOLS,
            tooltips=tooltips,
            toolbar_location=toolbar_location,
            title=title,
            max_width=width, height=height, 
        )
        p.toolbar.autohide = True
        return p
    
    def _get_html(self, p):
        html = bokeh.embed.file_html(p, bokeh.resources.CDN, "")
        return html


class DrugScreenDashboard(DataDashboard):
    
    def __init__(self, screen, treated, untreated, t0='T0', threshold=3, ctrl_label='negative_control',run_name='auto'):
        self.screen = screen
        self.threshold = threshold
        self.ctrl_label = ctrl_label
        self.run_name = run_name
        self.gamma_score_name = f'gamma:{untreated}_vs_{t0}'
        self.rho_score_name = f'rho:{treated}_vs_{untreated}'
        self.plots = {}
        super().__init__()

    def _prep_data(self,screen):

        gamma = screen.getPhenotypeScores(
            run_name=self.run_name,
            score_name=self.gamma_score_name,
            ctrl_label=self.ctrl_label,
            threshold=self.threshold,
        )

        rho = screen.getPhenotypeScores(
            run_name=self.run_name,
            score_name=self.rho_score_name,
            ctrl_label=self.ctrl_label,
            threshold=self.threshold
        )

        df = pd.DataFrame({
            'target': rho['target'],
            'rho_score': rho['score'],
            'rho_pvalue': rho['pvalue'],
            'rho_label': rho['label'],
            '-log10(rho_pvalue)': np.log10(rho['pvalue']) * -1,
            'gamma_score': gamma.loc[rho.index,'score'],
            'gamma_pvalue': gamma.loc[rho.index,'pvalue'],
            'gamma_label': gamma.loc[rho.index,'label'],
            '-log10(gamma_pvalue)': np.log10(gamma.loc[rho.index,'pvalue']) * -1,
        })

        return df
    
    def _plot_scatter(
            self,
            x_source,y_source,
            xaxis_label,yaxis_label,
            up_hit, down_hit,
            hit_label_col,
            x_min, x_max, y_min, y_max,
            title='',
            dot_size=1,
            width=500, height=400,
            toolbar_location='below',
            legend_loc="top_left"
        ):

        df = self._prep_data(self.screen)

        if y_max == 'auto': y_max = df[y_source].max() * 1.2
        if x_max == 'auto': x_max = df[x_source].max() * 1.2
        if y_min == 'auto': y_min = df[y_source].min() * 1.2
        if x_min == 'auto': x_min = df[x_source].min() * 1.2

        TOOLTIPS = [
            ("name", "@target"),
            ("rho score", "@rho_score"),
            ("rho p-value", "@rho_pvalue"),
            ("rho label", "@rho_label"),
            ("gamma score", "@gamma_score"),
            ("gamma p-value", "@gamma_pvalue"),
            ("gamma label", "@gamma_label"),        
        ]

        p = self._new_plot(
            title=title,
            tooltips=TOOLTIPS,
            width=width,
            height=height,
            toolbar_location=toolbar_location
        )

        source = bokeh.models.ColumnDataSource(
            df.loc[df[hit_label_col] == 'target_non_hit',:]
        )
        p.scatter(
            x=x_source, y=y_source,
            source=source,
            alpha=0.2,
            size=dot_size * 1.2, 
            color='gray', 
            legend_label='target_non_hit',
            name='circles'
        )

        # size_mapper=bokeh.models.LinearInterpolator(
        #     x=[df['1/gamma_score'].min(),df['1/gamma_score'].max()],
        #     y=[1,100]
        # )

        source = bokeh.models.ColumnDataSource(
            df.loc[df[hit_label_col] == up_hit,:]
        )
        p.scatter(
            x=x_source, y=y_source,
            source=source,
            alpha=0.8,
            size=dot_size * 1.2, 
            # size={'field':'1/gamma_score','transform':size_mapper},
            color='#fcae91',
            legend_label=up_hit,
            name='circles'
        )
        
        source = bokeh.models.ColumnDataSource(
            df.loc[df[hit_label_col] == down_hit,:]
        )
        p.scatter(
            x=x_source, y=y_source,
            source=source,
            alpha=0.8,
            # size={'field':'1/gamma_score','transform':size_mapper},
            size=dot_size * 1.2, 
            color='#bdd7e7',
            legend_label=down_hit,
            name='circles'
        )
        
        source = bokeh.models.ColumnDataSource(
            df.loc[df[hit_label_col] == self.ctrl_label,:]
        )
        p.scatter(
            x=x_source, y=y_source,
            source=source,
            alpha=0.2,
            size=dot_size*0.8, 
            color='silver',
            legend_label=self.ctrl_label,
            name='circles'
        )
        
        # Set x-axis and y-axis labels
        p.xaxis.axis_label = xaxis_label
        p.xaxis.axis_label_text_font_style = 'normal'
        p.yaxis.axis_label = yaxis_label
        p.yaxis.axis_label_text_font_style = 'normal'

        # Set x-axis limits
        p.x_range.start = x_min
        p.x_range.end = x_max

        # Set y-axis limits
        p.y_range.start = y_min
        p.y_range.end = y_max

        # Add legend
        if legend_loc == False or legend_loc == None:
            p.legend.visible = False
        else:
            p.legend.location = legend_loc

        p.title.text = title
        p.title.align = 'center'
        p.title.text_font_size = '12pt'
        p.title.text_font_style = 'bold'

        return p
    
    def RhoVolcanoPlot(
        self,
        x_source='rho_score', y_source='-log10(rho_pvalue)',
        xaxis_label='phenotype score',
        yaxis_label='-log10(p-value)',
        up_hit='resistance_hit', down_hit='sensitivity_hit',
        hit_label_col='rho_label',
        x_min=-2.5, x_max=2.5, y_min=0, y_max='auto',
        return_html=True,
        **kwargs
        ):
        p = self._plot_scatter(
            x_source, y_source,
            xaxis_label, yaxis_label,
            up_hit, down_hit,
            hit_label_col,
            x_min, x_max, y_min, y_max,
            **kwargs
        )

        if return_html:
            return self._get_html(p)
                    
        self.plots.update(
            {'RhoVolcanoPlot': p}
        )

    def GammaVolcanoPlot(
        self,
        x_source='gamma_score', y_source='-log10(gamma_pvalue)',
        xaxis_label='phenotype score',
        yaxis_label='-log10(p-value)',
        up_hit='up_hit', down_hit='essential_hit',
        hit_label_col='gamma_label',
        x_min=-2.5, x_max=2.5, y_min=0, y_max='auto',
        return_html=True,
        **kwargs
        ):
        p = self._plot_scatter(
            x_source, y_source,
            xaxis_label, yaxis_label,
            up_hit, down_hit,
            hit_label_col,
            x_min, x_max, y_min, y_max,
            **kwargs
        )

        if return_html:
            return self._get_html(p)
        
        self.plots.update(
            {'GammaVolcanoPlot': p}
        )
    
    def RhoGammaScatter(
        self,
        x_source='rho_score', y_source='gamma_score',
        xaxis_label='rho score',
        yaxis_label='gamma score',
        up_hit='resistance_hit', down_hit='sensitivity_hit',
        hit_label_col='rho_label',
        return_html=True,
        x_min=-2.5, x_max=2.5, y_min=-2.5, y_max=2.5,
        **kwargs
        ):
        p = self._plot_scatter(
            x_source, y_source,
            xaxis_label, yaxis_label,
            up_hit, down_hit,
            hit_label_col,
            x_min, x_max, y_min, y_max,
            **kwargs
        )

        if return_html:
            return self._get_html(p)
        
        self.plots.update(
            {'GammaRhoScatter': p}
        )
