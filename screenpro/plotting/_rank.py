import pandas as pd
import matplotlib.pyplot as plt
from ._utils import yellow_blue


def rank_plot(df, rank_col, color_col=None, name_col='target', highlight_values_dict=None, xlabel='Rank', ylabel='Values', title='Rank Plot', ax=None, dot_size=1.5, highlight_size_factor=100, txt_font_size=8, **args):
    """
    Plot the ranks against their values with specified color.

    Args:
        df (DataFrame): The input DataFrame.
        rank_col (str): The column name containing the values to be ranked.
        color_col (str): The column name containing the values to be used for color coding. Default is None.
        name_col (str, optional): The column name containing the names of the values. Default is 'target'.
        highlight_values_dict (dict, optional): A dictionary specifying the values to be highlighted. 
            The keys are the highlight colors and the values are dictionaries with 'genes' and 'text' keys. 
            'genes' is a list of values to be highlighted and 'text' is a boolean indicating whether to display 
            the names of the highlighted values. Default is None.
        xlabel (str, optional): The label for the x-axis. Default is 'Rank'.
        ylabel (str, optional): The label for the y-axis. Default is 'Values'.
        title (str, optional): The title of the plot. Default is 'Rank Plot'.
        ax (matplotlib.axes.Axes, optional): The axis object to plot on. If not provided, a new axis will be created.
        dot_size (float, optional): The size of the dots in the scatter plot. Default is 1.5.
        highlight_size_factor (int, optional): The size factor for the highlighted dots. Default is 100.
        txt_font_size (int, optional): The font size for the text labels. Default is 8.
        **args: Additional keyword arguments to be passed to the scatter plot.

    Returns:
        matplotlib.axes.Axes: The axis object containing the plot.
    """
    # Create a new DataFrame with the values and their corresponding ranks
    rank_df = df.copy()
    rank_df['Rank'] = rank_df[rank_col].rank()
    rank_df.sort_values('Rank', inplace=True)

    # Use a color that is suitable for publications
    if color_col is None:
        color_col = 'darkgray'

    # If no axis is provided, create one
    if ax is None:
        _, ax = plt.subplots()

    # Plot the ranks against their values with specified color
    rank_df.plot.scatter(
        'Rank', rank_col, marker='o',
        colormap=yellow_blue,
        s=dot_size,
        c=color_col, ax=ax,
        colorbar=False,
        **args
    )

    if highlight_values_dict is not None:
        for highlight_color, highlight_values in highlight_values_dict.items():
            highlight_ranks = rank_df[rank_df[name_col].isin(highlight_values['genes'])]
            ax.plot(highlight_ranks['Rank'], highlight_ranks[rank_col], 'o', color=highlight_color, markersize=dot_size * highlight_size_factor)
    
            if highlight_values['text'] is not False:
                for i, row in highlight_ranks.iterrows():
                    ax.text(row['Rank'] + .01, row[rank_col] + .001, row[name_col], fontsize=txt_font_size, color=highlight_color, ha='right')

    # Add labels and title
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    # Customize the grid lines for a clean look
    ax.grid(False)

    return rank_df, ax