#!/usr/bin/env python3


"""
Functions to be used to generate plots.
"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pathlib as Path
import seaborn as sns


def main():

    return None


def fastqc_plots(
    fastqc_path: Path,
    output_path: Path,
    fqc_str: str = 'fastqc',
    fqc_ext: str = '.zip'
) -> None:
    """
    Handles the generation of all plots associated with FastQC metrics.
    """
    # Get all files within dir_path
    for file in fastqc_path.iterdir():
        if file.is_file() and \
            fqc_str.lower() in file.name.lower() and \
            file.name.endswith(fqc_ext):
            pass
    return None


def _plot_boxplot(
    data_bp: list[np.ndarray],
    output_path: Path = '',
    show_plot: bool = True
) -> None:
    """
    Generates a boxplot from a list of groups of data points.

    Args:
        data_df (pd.DataFrame): DataFrame of data points to be plotted.
        output_path (Path): Path to the output file where the plot is 
                            saved.
        show_plot (bool): Defines if the plot is shown or not.
    """
    # Create the boxplot
    plt.boxplot(data_bp)

    # Add a title and labels
    plt.title("Boxplot")
    plt.xlabel("Data Sets")
    plt.ylabel("Values")

    # Save the plot if output_path is defined
    if output_path:
        plt.savefig(output_path)

    # Display the plot if show_plot is True
    if show_plot:
        plt.show()
    
    return None


def _plot_heatmap(
    data_df: pd.DataFrame,
    output_path: Path = '',
    show_plot: bool = True
) -> None:
    """
    Generates a heatmap from a grid of data points.

    Args:
        data_df (pd.DataFrame): DataFrame of data points to be plotted.
        output_path (Path): Path to the output file where the plot is 
                            saved.
        show_plot (bool): Defines if the plot is shown or not.
    """

    plt.figure(figsize=(8, 6)) # Set the figure size
    sns.heatmap(data_df, annot=True, cmap='coolwarm', fmt=".2f", linewidths=.5, vmin=-1, vmax=1)
    plt.title('Heatmap')

    # Save the plot if output_path is defined
    if output_path:
        plt.savefig(output_path)

    # Display the plot if show_plot is True
    if show_plot:
        plt.show()
    
    return None


def _plot_histogram(
    data_series: pd.Series,
    output_path: Path = '',
    show_plot: bool = True
) -> None:
    """
    Generates a histogram from a series of data points.

    Args:
        data_series (pd.Series): Series of data points to be plotted.
        output_path (Path): Path to the output file where the plot is 
                            saved.
        show_plot (bool): Defines if the plot is shown or not.
    """

    # Create the histogram
    plt.hist(data_series, bins=30, color='skyblue', edgecolor='black')

    # Add labels and a title for clarity
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title('Histogram')

    # Save the plot if output_path is defined
    if output_path:
        plt.savefig(output_path)

    # Display the plot if show_plot is True
    if show_plot:
        plt.show()
    
    return None


if __name__ == "__main__":
    main()
