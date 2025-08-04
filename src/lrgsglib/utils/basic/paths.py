from .common import *
from .iterables import uniques
#
__all__ = [
    'list_dir',
    'remove_directory_if_empty',
    'remove_empty_dirs',
    'find_matching_files',
    'extract_value_from_filename',
    'extract_values_from_filenames',
    'extract_and_sort_values',
]
#
logger = logging.getLogger(__name__)
#
def list_dir(path: str) -> List[str]:
    """
    Return a sorted list of names of all entries in the specified directory.

    Parameters
    ----------
    path : str
        Path to the target directory.

    Returns
    -------
    List[str]
        Sorted list of entry names (files and subdirectories).

    Raises
    ------
    FileNotFoundError
        If the specified path does not exist.
    NotADirectoryError
        If the specified path is not a directory.
    PermissionError
        If access to the directory is denied.

    Example
    -------
    >>> list_dir("/home/user/docs")
    ['file1.txt', 'file2.txt', 'reports']
    """
    dir_path = Path(path)
    if not dir_path.exists():
        raise FileNotFoundError(f"Path not found: {path}")
    if not dir_path.is_dir():
        raise NotADirectoryError(f"Not a directory: {path}")
    return sorted(entry.name for entry in dir_path.iterdir())
#
def remove_directory_if_empty(path: Union[str, Path]) -> bool:
    """
    Remove a directory (and all its empty subdirectories) if it contains no files.

    Parameters
    ----------
    path : Union[str, Path]
        Path to the directory to check.

    Returns
    -------
    bool
        True if the directory existed and was removed; False otherwise
        (not a directory, contains files, or removal failed).
    """
    p = Path(path)
    if not p.is_dir():
        return False

    # If any file exists anywhere under p, do not remove
    if any(item.is_file() for item in p.rglob('*')):
        return False

    # Safe to remove entire tree of empty dirs
    try:
        shutil.rmtree(p)
        return True
    except OSError as e:
        logger.error(f"Failed to remove {p!r}: {e}")
        return False
#
def remove_empty_dirs(path: Union[str, Path]) -> None:
    """
    Recursively remove all empty subdirectories under the given directory.

    Walks the directory tree bottom-up, attempting to remove each subdirectory
    that is empty. Errors (e.g. non-empty or permission issues) are logged
    but do not interrupt the process.

    Parameters
    ----------
    path : Union[str, Path]
        Root directory under which to remove empty subdirectories.

    Returns
    -------
    None
    """
    root = Path(path)
    if not root.is_dir():
        return

    # Sort by depth (deepest first) so children are removed before parents
    dirs = [p for p in root.rglob('*') if p.is_dir()]
    dirs.sort(key=lambda p: len(p.parts), reverse=True)

    for d in dirs:
        try:
            d.rmdir()
            logger.info(f"Removed empty directory: {d}")
        except OSError:
            # Not empty or cannot remove; skip
            continue

def find_matching_files(search_dir: str, pattern_str: str) -> Union[str, List[str]]:
    """
    Searches all files in the specified directory that match a given pattern.

    Parameters
    ----------
    search_dir : str
        The directory in which to search for files.
    pattern_str : str
        The string pattern to search for within the file names. This pattern
        will be compiled into a regular expression.

    Returns
    -------
    Union[str, List[str]]
        A single file name if only one file matches the pattern, or a list of 
        file names within the search directory that match the given pattern.
    """
    pattern = re.compile(f'.*{pattern_str}.*')
    all_files = list_dir(search_dir)
    matching_files = [fname for fname in all_files if pattern.match(fname)]
    
    if len(matching_files) == 1:
        return matching_files[0]
    return matching_files


def extract_value_from_filename(file_name: str, value_pattern: str) -> float:
    """
    Extracts a numeric value from a single file name based on a specified regular expression pattern.
    
    Parameters
    ----------
    file_name : str
        The file name from which to extract the value.
    value_pattern : str
        The regular expression pattern used to extract the numeric value from the file name. The pattern
        should contain a capturing group for the numeric value.
        
    Returns
    -------
    float
        The extracted numeric value from the file name.
        
    Examples
    --------
    >>> file_name = "data_p=1.5.pkl"
    >>> value_pattern = r"p=([\\d.]+)"
    >>> extract_value_from_filename(file_name, value_pattern)
    1.5
    
    This function extracts the p-value from the provided file name.
    """
    match = re.search(value_pattern, file_name)
    if match:
        return float(match.group(1).rstrip('.'))
    else:
        raise ValueError("No match found in the file name.")


def extract_values_from_filenames(file_names: List[str], value_pattern: str, sort: bool = True, unique: bool = False) -> NDArray:
    """
    Extracts numeric values from a list of file names based on a specified regular expression pattern.
    
    Parameters
    ----------
    file_names : List[str]
        A list of file names from which to extract the values.
    value_pattern : str
        The regular expression pattern used to extract numeric values from the file names. The pattern
        should contain a capturing group for the numeric value.
    sort : bool, optional
        Specifies whether to sort the extracted values. Default is True.
        
    Returns
    -------
    numpy.ndarray
        An array of extracted numeric values from the file names. If `sort` is True, this array will be
        sorted in ascending order.
        
    Examples
    --------
    >>> file_names = ["data_p=1.5.pkl", "experiment_p=2.0.pkl", "results_p=0.5.pkl"]
    >>> value_pattern = "p=([\\d.]+)"
    >>> extract_values_from_filenames(file_names, value_pattern)
    array([0.5, 1.5, 2.0])
    
    This function extracts the p-values from the provided file names and returns them sorted in ascending order.
    """
    # values = [re.search(value_pattern, filename).group(1) 
    #           for filename in file_names if re.search(value_pattern, filename)]
    values = [
        float(match.group(1).rstrip('.'))
        for file_name in file_names
        if (match := re.search(value_pattern, file_name))
    ]
    if unique:
        values = uniques(values)
    if sort:
        values.sort()
    return np.array(values)

def extract_and_sort_values(path: str, search_pattern: str, value_pattern: str = None, sort: bool = True) -> NDArray:
    """
    Searches for files in a given directory that match a specified pattern, extracts numerical values from
    these file names based on another pattern (assumed to be the same as search pattern if not provided), 
    and optionally sorts these values. Raises an error if value_pattern does not contain a capturing group.
    
    Parameters
    ----------
    path : str
        The directory path to search for files.
    search_pattern : str
        The regular expression pattern to match file names for the search.
    value_pattern : str, optional
        The regular expression pattern used to extract numeric values from the matched file names. If None, 
        search_pattern is used. This pattern must contain at least one capturing group.
    sort : bool, optional
        Specifies whether to sort the extracted values. Default is True.
        
    Returns
    -------
    numpy.ndarray
        An array of extracted numeric values from the matched file names. If `sort` is True, this array will be
        sorted in ascending order.
        
    Raises
    ------
    ValueError
        If value_pattern does not contain at least one capturing group.
    """
    if value_pattern is None:
        value_pattern = search_pattern

    # Check if the value_pattern contains at least one capturing group
    if not re.search(r"\((?!\?:).+?\)", value_pattern):
        raise ValueError("value_pattern must contain at least one capturing group.")
    
    file_names = find_matching_files(path, search_pattern)
    return extract_values_from_filenames(file_names, value_pattern, sort)