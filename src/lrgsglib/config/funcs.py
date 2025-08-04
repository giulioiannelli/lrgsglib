#
import glob
import re
#
import numpy as np
#
from numpy.typing import NDArray
#
from collections import Counter, defaultdict
from os import chdir
from pathlib import Path
from typing import Optional, Union, List, Tuple, Dict
#
from ..config.const import *

__all__ = [
    'move_to_rootf',
    'peq_fstr',
    'Teq_fstr',
    'avgeq_fstr',
    'build_pT_fname',
    'build_fname_or_pattern_direct',
    'read_files_to_2d_array',
    'find_shared_p_values',
    'bin_eigenvalues',
]
#
def move_to_rootf(print_tf: bool = True, pathname: str = None):
    """
    Move to the root directory of the current working directory.

    Parameters:
    -----------
    print_tf : bool, optional
        If True, print the current working directory after moving to the root.
        Default is True.

    Notes:
    ------
    - The function continuously moves up the directory hierarchy ('../') until it reaches
      a directory with the name specified by the 'PATH_ROOTF' constant.
    - If 'print_tf' is set to True, it prints the current working directory after the move.

    Example:
    --------
    To move to the root directory of the current working directory and print the path:
    move_to_rootf(print_tf=True)
    """
    pcwd = Path.cwd()
    if pathname is None:
        pathname = PATHNLLIB
    try:
        while Path.cwd().name != pathname:
            chdir('../')
            if Path.cwd().name == '/':
                break
        if Path.cwd().name == '/':
            raise FileNotFoundError(f"Root directory '{pathname}' not found.")
    except FileNotFoundError as e:
        chdir(pcwd)
        print(e)
    if print_tf:
        print("Current working directory:", Path.cwd())


def peq_fstr(p: float) -> str:
    """
    Format the probability of flip (p) as a string.

    Parameters
    ----------
    p : float
        Probability of flip.

    Returns
    -------
    str
        Formatted string representing the probability of flip.

    Example
    -------
    >>> peq_fstr(0.25)
    'p=0.25'
    """
    return f"p={p:.3g}"

def Teq_fstr(T: float) -> str:
    """
    Format the temperature (T) as a string.

    Parameters
    ----------
    T : float
        Temperature value.

    Returns
    -------
    str
        Formatted string representing the temperature.

    Example
    -------
    >>> Teq_fstr(1.5)
    'T=1.5'
    """
    return f"T={T:.3g}"

def avgeq_fstr(avg: int) -> str:
    """
    Format the number of averages as a string.

    Parameters
    ----------
    avg : int
        Number of averages.

    Returns
    -------
    str
        Formatted string representing the number of averages.

    Example
    -------
    >>> avgeq_fstr(10)
    'avg=10'
    """
    return f"avg={avg:d}"

def build_p_fname(
    base: str,
    pflip: float,
    out_suffix: str = '',
    path: Optional[Union[str, Path]] = None,
    ext: str = '.bin'
) -> Union[str, Path]:
    from ..utils.basic.strings import join_non_empty
    fname = join_non_empty('_', base, peq_fstr(pflip), out_suffix) + ext
    return Path(path) / fname if path else Path(fname)

def build_pT_fname(
    base: str,
    pflip: float,
    T: float,
    out_suffix: str = '',
    path: Optional[Union[str, Path]] = None,
    ext: str = ''
) -> Union[str, Path]:
    """
    Build a filename for storing or loading data, based on parameters.

    Args:
        base (str): Base name for the file.
        pflip (float): Probability of flip, used in filename.
        T (float): Temperature, used in filename.
        out_suffix (str, optional): Output suffix to append to filename. Defaults to ''.
        path (str or Path, optional): Directory path to prepend to filename. Defaults to None.
        ext (str, optional): File extension. Defaults to '.bin'.

    Returns:
        str or Path: The constructed filename, as a Path if path is provided, else as a string.
    """
    from ..utils.basic.strings import join_non_empty
    fname = join_non_empty('_', base, peq_fstr(pflip), Teq_fstr(T), out_suffix) + ext
    if path:
        return Path(path) / fname
    else:
        return fname

def build_fname_or_pattern_direct(
    basename: str,
    p: float,
    T: float,
    out_suffix: str,
    mode: str = 'fname',
    n_avg: int = 0,
    ext: str = ''
) -> str:
    """
    Build a filename or pattern for saving results based on direct parameters.

    Parameters
    ----------
    basename : str
        The base name for the file, typically including the operation mode (e.g., 'smatch_recon').
    p : float
        The 'p' parameter value.
    T : float
        The 'T' parameter value.
    out_suffix : str
        An output suffix to append to the filename base.
    mode : str, optional
        Mode for the filename, either 'fname' or 'pattern', by default 'fname'.
    n_avg : int, optional
        Number of averages to include in the filename if mode is 'fname', by default 0.
    ext : str, optional
        Extension to be used by build_pT_fname for the fname_base part, by default ''.

    Returns
    -------
    str
        The constructed filename or pattern.
    """
    fname_base = build_pT_fname(basename, p, T, out_suffix=out_suffix)
    match mode:
        case 'fname':
            # Ensure n_avg is handled correctly, avgeq_fstr might expect non-zero for meaningful output
            avg_str = avgeq_fstr(n_avg) if n_avg > 0 else ''
            from ..utils.basic.strings import join_non_empty
            fname = join_non_empty('_', fname_base, avg_str)
        case 'pattern':
            fname = '_'.join([fname_base, 'avg=*'])
        case _:
            raise ValueError(f"Unsupported mode '{mode}' for filename construction.")
    return fname + ext if ext else fname















def read_files_to_2d_array(folder_path, keyword):
    """
    Read files from a folder that contain a specific keyword in their name and append each file's contents to a 2D array.

    Parameters
    ----------
    folder_path : str
        Path to the folder containing the files.
    keyword : str
        String that must be part of the file's name to be processed.

    Returns
    -------
    numpy.ndarray
        A 2D array containing the contents of each processed file, with the data type of the array elements as float.
    """
    # Initialize the 2D array
    data_2d_array = []
    
    # List all files in the given folder
    from ..utils.basic.paths import list_dir
    for file_name in list_dir(folder_path):
        # Check if the file name contains the keyword
        if keyword in file_name.split('_'):
            # Construct full file path
            file_path = Path(folder_path) / file_name
            # Open and read the file
            with open(file_path, 'r') as file:
                # Assuming each line of a file represents a separate data entry
                file_contents = [line.strip() for line in file.readlines()]
                data_2d_array.append(file_contents)
    
    return np.array(data_2d_array).astype(float)


def find_shared_p_values(pattern: str, pathdir: str, extension: str = '.pkl') -> List[Tuple[float, int]]:
    """
    Identifies and counts p-values found in file paths that appear in at least two different subdirectories,
    considering files with a specific extension.
    
    This function searches for files with the given extension within subdirectories of a specified path, 
    extracts p-values based on a provided pattern, and counts the occurrences across different subdirectories. 
    It returns a sorted list of unique p-values that are present in two or more subdirectories, along with 
    their occurrence count.

    Parameters
    ----------
    pattern : str
        The regular expression pattern used to extract p-values from file paths. 
        The pattern should contain a capturing group for the p-value.
    pathdir : str
        The root directory path under which to search for files with the specified extension. 
        This path is expected to contain subdirectories where the files are located.
    extension : str, optional
        The file extension to look for (default is '.pkl').

    Returns
    -------
    list of tuple
        A sorted list of tuples, where each tuple contains a p-value and the count of subdirectories 
        in which the p-value is found. The list is sorted by p-value in ascending order.
        
    Examples
    --------
    >>> find_shared_p_values(r"p=([\\d.]+)", "/path/to/data/")
    [(1.0, 2), (2.5, 3)]
    
    This example indicates that p-value 1.0 was found in 2 subdirectories, and p-value 2.5 was found in 3 subdirectories.
    """
    
    # Dictionary to hold sets of subdirectories for each found p value
    p_values_dirs = defaultdict(set)

    # Use glob to iterate over all files with the specified extension in subfolders
    for filepath in glob.glob(f'{pathdir}*/*{extension}'):
        match = re.search(pattern, filepath)
        if match:
            # Extract the p value
            p_value = float(match.group(1))
            # Extract subdirectory from the filepath
            # Adjust the split index based on your path structure
            subdirectory = filepath.split('/')[4]
            # Add the subdirectory to the set for this p value
            p_values_dirs[p_value].add(subdirectory)

    # Prepare a list to hold p values and the number of sharing subdirectories
    p_values_shared_count = []

    # Filter and count p values that appear in at least two different subdirectories
    for p_value, dirs in p_values_dirs.items():
        num_shared = len(dirs)
        if num_shared >= 2:
            p_values_shared_count.append((p_value, num_shared))

    # Sort the list by p value
    p_values_shared_count.sort()
    
    return p_values_shared_count







def bin_eigenvalues(eigenvalues: List[float], bins: np.ndarray, bin_centers: np.ndarray) -> Counter:
    """
    Bins eigenvalues into specified bins and returns the count for each bin using bin centers as keys.

    Parameters
    ----------
    eigenvalues : List[float]
        The eigenvalues to be binned.
    bins : np.ndarray
        The edges of the bins.
    bin_centers : np.ndarray
        The center values of each bin.

    Returns
    -------
    Counter
        A dictionary-like collection where keys are bin centers and values are the count of eigenvalues in each bin.
    """
    bin_indices = np.digitize(eigenvalues, bins, right=True) - 1
    bin_indices = np.clip(bin_indices, 0, len(bin_centers) - 1)  # Ensure indices are within the valid range
    bin_keys = [bin_centers[index] for index in bin_indices]
    return Counter(bin_keys)












