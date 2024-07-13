import pandas as pd
import matplotlib.pyplot as plt


def rankPlot(series, ax=None, highlight_values_dict=None, xlabel='Rank', ylabel='Values', title='Rank Plot', dot_size=3, highlight_size_factor=100):
    # Create a copy of the series and sort it
    sorted_series = series.sort_values(ascending=True)

    # Create a new DataFrame with the values and their corresponding ranks
    rank_df = pd.DataFrame({'Rank': range(1, len(sorted_series) + 1), 'Values': sorted_series.values, 'Name': sorted_series.index})

    # Use a color that is suitable for publications
    plot_color = 'darkgray'

    # If no axis is provided, create one
    if ax is None:
        fig, ax = plt.subplots()

    # Plot the ranks against their values with specified color
    ax.plot(rank_df['Rank'], rank_df['Values'], marker='o', linestyle='-', color=plot_color, markersize=dot_size)

    if highlight_values_dict is not None:
        for highlight_color, highlight_values in highlight_values_dict.items():
            highlight_ranks = rank_df[rank_df['Name'].isin(highlight_values['genes'])]
            ax.plot(highlight_ranks['Rank'], highlight_ranks['Values'], 'o', color=highlight_color, markersize=dot_size * highlight_size_factor)
    
            if highlight_values['text'] is not False:
                for i, row in highlight_ranks.iterrows():
                    ax.text(row['Rank'] + .01, row['Values'] + .001, row['Name'], fontsize=8, color=highlight_color, ha='right')

    # Add labels and title
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    # Customize the grid lines for a clean look
    ax.grid(False)

    return ax
