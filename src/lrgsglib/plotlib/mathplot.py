from ..utils.basic import log_binning
from .const_plotlib import *
#
__all__ = [
    'plot_log_distribution',
]
#
def plot_log_distribution(data, fig_ax=None, binnum=20, 
                          xlabel="X-Axis (log scale)", ylabel="Frequency", 
                          title="Log-Scale Distribution", grid=True, 
                          log_scale=True, **kwargs):
    """
    Plots the distribution of the given data using a logarithmic scale on the x-axis.

    Parameters:
    -----------
    data : np.ndarray or list
        The data to be plotted in the distribution.

    fig_ax : tuple, optional
        Tuple of (fig, ax) to use an existing figure and axis. If None, a new 
        figure and axis are created. Default is None.

    binnum : int, optional
        Number of bins for the histogram. Default is 20.

    xlabel : str, optional
        Label for the x-axis. Default is "X-Axis (log scale)".

    ylabel : str, optional
        Label for the y-axis. Default is "Frequency".

    title : str, optional
        Title for the plot. Default is "Log-Scale Distribution".

    grid : bool, optional
        Whether to display a grid. Default is True.

    log_scale : bool, optional
        Whether to apply a logarithmic scale to the x-axis. Default is True.

    **kwargs : optional
        Additional keyword arguments for the `plt.bar` function to customize the
        plot appearance (e.g., `color`, `alpha`, etc.).

    Returns:
    --------
    None
        Displays the plot using the provided or created figure and axis.
    """
    # Logarithmic binning of data
    bin_centers, hist, bin_w = log_binning(data, binnum)
    
    # Check if custom fig and ax are provided, otherwise create them
    if fig_ax is None:
        fig, ax = plt.subplots()
    else:
        fig, ax = fig_ax
    
    # Plot the histogram using a bar chart with customizable kwargs
    ax.bar(bin_centers, hist, width=bin_w, alpha=0.75, ec='black', **kwargs)
    
    # Apply logarithmic scale if specified
    if log_scale:
        ax.set_xscale('log')
    
    # Set labels, title, and grid as per the provided arguments
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    if grid:
        ax.grid(True, which="both", ls="--")  # Grid for both major and minor ticks
    
    # Show the plot if no custom fig_ax is provided
    if fig_ax is None:
        plt.show()