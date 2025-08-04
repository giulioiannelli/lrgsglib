from .common import *
#
__all__ = [
    'boolean_overlap_fraction',
    'cProd_Iter',
    'cProdSel_Iter',
    'cProd_Iter_adj',
    'extract_subdictionary',
    'first_index_changing_condition',
    'flatten',
    'inf_array_regularization',
    'sort_array_by_column',
    'sum_tuples',
    'uniques',
    'unzip_dict_items'
]
#
def boolean_overlap_fraction(boolist1, boolist2):
    """
    Calculate the fraction of overlapping True values between two boolean lists.

    Parameters:
    -----------
    boolist1 : list of bool
        The first boolean list.

    boolist2 : list of bool
        The second boolean list.

    Returns:
    --------
    float
        The fraction of overlapping True values between the two boolean lists.

    Notes:
    ------
    - The function computes the overlap fraction by first performing a bitwise XOR (^)
      operation between the two boolean lists, which results in a new boolean list where
      True represents differences and False represents matches.
    - The bitwise NOT (~) operator is then applied to invert the differences,
      turning them into True values.
    - The 'sum' function counts the number of True values in the inverted list,
      representing the number of overlapping True values.
    - Finally, this count is divided by the length of 'boolist1' to obtain the overlap fraction.

    Example:
    --------
    boolist1 = [True, False, True, False, True]
    boolist2 = [True, True, False, False, True]
    overlap_fraction = boolean_overlap_fraction(boolist1, boolist2)
    # The result is 0.4, indicating 40% overlap of True values between the two lists.

    """
    return sum(~(boolist1 ^ boolist2))/len(boolist1)
#
def cProd_Iter(dim: Union[int, Tuple]) -> iter:
    """
    Generates the Cartesian product for an n-dimensional space.

    This function creates an iterable over all possible combinations of coordinates
    in an n-dimensional grid, defined by the dimensions specified in the input tuple.
    Each element in the 'dim' tuple represents the size of the grid along that dimension.
    The Cartesian product is generated using a sequence of range objects, one for each
    dimension, thereby producing a sequence of tuples that represent points in the
    n-dimensional space.

    Parameters
    ----------
    dim : tuple of int
        A tuple where each element specifies the size of the grid along a particular
        dimension. For example, (2, 3) represents a 2D grid with dimensions 2x3.

    Returns
    -------
    iterator
        An iterator over tuples, where each tuple represents a point in the
        n-dimensional space defined by 'dim'. Each tuple contains integers,
        with the i-th integer representing the coordinate along the i-th dimension.

    Examples
    --------
    >>> list(cartesian_product((2, 3)))
    [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]

    >>> list(cartesian_product((1, 2, 3)))
    [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 0), (0, 1, 1), (0, 1, 2)]

    Notes
    -----
    - The order of the tuples in the output iterator follows the lexicographical order
      based on the input dimensions. This means smaller coordinates come before larger ones.
    - This function is particularly useful for iterating over multi-dimensional arrays
      or grids, where you need to visit each cell or point in the space.
    - The implementation relies on itertools.product, which is efficient and avoids
      explicitly constructing the grid in memory, making it suitable for large dimensions.

    See Also
    --------
    itertools.product : Cartesian product of input iterables.
    """
    return product(*[range(d) for d in dim])
#
def cProdSel_Iter(dim: Union[int, Tuple], selected_indices: Union[List, Tuple]) -> iter:
    """
    Generates the Cartesian product for selected dimensions in an n-dimensional space.

    This function creates an iterable over all possible combinations of coordinates
    in the specified dimensions of an n-dimensional grid. It allows for focusing on
    a subset of all available dimensions, which is specified by the selected_indices.

    Parameters
    ----------
    dim : tuple of int
        A tuple where each element specifies the size of the grid along a particular
        dimension. For example, (2, 3, 4) represents a 3D grid with dimensions 2x3x4.
    selected_indices : list or tuple of int
        Indices of the dimensions to include in the Cartesian product. For example,
        (1, 2) selects the second and third dimensions from 'dim'.

    Returns
    -------
    iterator
        An iterator over tuples, where each tuple represents a point in the
        specified dimensions of the n-dimensional space. Each tuple contains integers,
        with the i-th integer representing the coordinate along the selected i-th dimension.

    Examples
    --------
    >>> list(cartesian_product_selected((2, 3, 4), (0, 2)))
    [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (1, 3)]

    >>> list(cartesian_product_selected((2, 3, 4, 5), (1, 3)))
    [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (2, 0), (2, 1), (2, 2), (2, 3), (2, 4)]

    Notes
    -----
    - The function is useful for iterating over a selected subset of dimensions in a multi-dimensional space.
    - By allowing specification of dimensions of interest, it provides flexibility for applications that do not require
      the full Cartesian product of all dimensions.
    """
    # Generate ranges for selected dimensions
    ranges = [range(dim[i]) for i in selected_indices]
    
    # Return Cartesian product of selected dimensions
    return product(*ranges)

def cProd_Iter_adj(dim: Union[int, Tuple], range_adjustment: Union[int, List] = 0) -> iter:
    """
    Generates the Cartesian product for an n-dimensional space with adjustable ranges.

    This function creates an iterable over all possible combinations of coordinates
    in an n-dimensional grid, defined by the dimensions specified in the input tuple,
    adjusted by the range_adjustment which can be a single integer or a list of integers.

    Parameters
    ----------
    dim : tuple of int
        A tuple where each element specifies the size of the grid along a particular
        dimension. For example, (2, 3) represents a 2D grid with dimensions 2x3.
    range_adjustment : int or list of int, optional
        An integer to adjust all dimensions equally, or a list of integers to adjust
        each dimension individually. The default value is 0, which means no adjustment.

    Returns
    -------
    iterator
        An iterator over tuples, where each tuple represents a point in the
        adjusted n-dimensional space defined by 'dim' and 'range_adjustment'.
        Each tuple contains integers, with the i-th integer representing the
        coordinate along the i-th dimension.

    Examples
    --------
    >>> list(cProd_Iter((2, 3), 1))
    [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (1, 3), (2, 0), (2, 1), (2, 2), (2, 3)]

    >>> list(cProd_Iter((1, 2), [1, 2]))
    [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (1, 3)]

    Notes
    -----
    The `range_adjustment` allows for flexible configuration of the generated space,
    accommodating scenarios that require padding or extending the dimensions.
    """
    # Handle both int and list types for range_adjustment
    if isinstance(range_adjustment, int):
        adjusted_ranges = [range(d + range_adjustment) for d in dim]
    elif isinstance(range_adjustment, list):
        adjusted_ranges = [range(d + adj) for d, adj in zip(dim, range_adjustment)]
    else:
        raise ValueError("range_adjustment must be an int or a list of ints")

    return product(*adjusted_ranges)
#
def extract_subdictionary(dictionary: dict, subkeys: list) -> dict:
    """
    Extract a subdictionary containing only the specified keys from the input dictionary.

    Parameters
    ----------
    dictionary : dict
        The input dictionary from which to extract key-value pairs.
    subkeys : list
        A list of keys to extract from the input dictionary.

    Returns
    -------
    dict
        A new dictionary containing only the key-value pairs where the keys are in `subkeys` and exist in the input dictionary.

    Notes
    -----
    - If a key in `subkeys` does not exist in the input dictionary, it is ignored.
    - This function is useful for filtering dictionaries to include only relevant keys.

    Examples
    --------
    >>> input_dict = {'a': 1, 'b': 2, 'c': 3}
    >>> keys_to_extract = ['a', 'c']
    >>> extract_subdictionary(input_dict, keys_to_extract)
    {'a': 1, 'c': 3}

    >>> extract_subdictionary(input_dict, ['d'])
    {}
    """
    # Use dictionary comprehension to filter the input dictionary
    return {key: dictionary[key] for key in subkeys if key in dictionary}
#
def first_index_changing_condition(condition):
    """
    Find the index of the first change in a boolean condition.

    Parameters:
    -----------
    condition : numpy.ndarray of bool
        A boolean condition represented as a NumPy array.

    Returns:
    --------
    int
        The index of the first change in the condition.

    Notes:
    ------
    - The function compares adjacent elements in the 'condition' array and returns
      the index of the first change from True to False or vice versa.
    - If there are no changes in the condition, the function returns 0.

    Example:
    --------
    condition = np.array([True, True, True, False, False, True])
    first_change_index = first_index_changing_condition(condition)
    # The result is 3, indicating the first change from True to False occurred at index 3.

    """
    return np.where(condition[:-1] != condition[1:])[0][0]
#
def flatten(xs):
    """
    Recursively flattens a nested iterable into a flat iterable.

    Parameters:
    -----------
    xs : iterable
        The input nested iterable to be flattened.

    Returns:
    --------
    object
        The flattened elements from the input iterable.

    Notes:
    ------
    This function takes a nested iterable and yields each element from the
    nested structure as a flat iterable. It recursively processes nested
    iterables such as lists or tuples.

    Example:
    --------
    >>> nested_list = [1, [2, [3, 4], 5], 6]
    >>> flattened = list(flatten(nested_list))
    >>> flattened
    [1, 2, 3, 4, 5, 6]
    """
    for x in xs:
        if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
            yield from flatten(x)
        else:
            yield x
#
def inf_array_regularization(arrinfs: NDArray) -> NDArray:
    """
    Regularizes an array by replacing infinite values with extremes of non-infinite values.

    Parameters:
    - arrinfs (ndarray): Input array containing infinite and finite numerical values.

    Returns:
    - ndarray: Regularized array with infinite values replaced by extremes of non-infinite values.

    Examples:
    >>> arr = np.array([1, 2, np.inf, 5, -np.inf, 7])
    >>> inf_array_regularization(arr)
    array([1., 2., 7., 5., 7., 7.])
    """
    # Filtering out infinite values from the input array (arrinfs)
    arrinfs_nnans = arrinfs[(arrinfs != np.inf) & (arrinfs != -np.inf)]

    # Regularizing infinite values in the array using nan_to_num function
    arrinfs = np.nan_to_num(arrinfs, posinf=np.max(arrinfs_nnans),
                            neginf=np.min(arrinfs_nnans))

    return arrinfs
#
def sort_array_by_column(arr: NDArray, column_index: int) -> NDArray:
    """
    Sorts a numpy array by a specific column.

    Parameters
    ----------
    arr : np.ndarray
        The array to be sorted.
    column_index : int
        The index of the column to sort by.

    Returns
    -------
    np.ndarray
        The sorted array.
    """
    return arr[arr[:, column_index].argsort()]
#
def sum_tuples(tuple1: tuple, tuple2: tuple) -> tuple:
    """
    Sum two tuples element-wise.

    If the tuples are of different lengths, the function sums elements up to the
    length of the shorter tuple, ignoring extra elements in the longer tuple.

    Parameters:
    -----------
    tuple1 : tuple
        The first tuple to be summed.
    tuple2 : tuple
        The second tuple to be summed.

    Returns:
    --------
    tuple
        A new tuple containing the element-wise sums of `tuple1` and `tuple2`.

    Example:
    --------
    >>> sum_tuples((1, 2, 3), (4, 5, 6))
    (5, 7, 9)
    """
    return tuple(a + b for a, b in zip(tuple1, tuple2))
#
def uniques(lst: List[Any]) -> List[Any]:
    """
    Returns a list of unique elements from the input list. This function 
    leverages Python's built-in `set` data structure to eliminate duplicate 
    entries efficiently. Note that the original order of elements is not 
    preserved.

    Parameters
    ----------
    lst : List[Any]
        The input list from which to extract unique elements. The list can 
        contain elements of any data type that is hashable.

    Returns
    -------
    List[Any]
        A new list containing only the unique elements from the input list, 
        with duplicates removed.

    Examples
    --------
    ```python
    >>> uniques([1, 2, 2, 3, 4, 4, 5])
    [1, 2, 3, 4, 5]

    >>> uniques(['apple', 'banana', 'apple', 'cherry'])
    ['apple', 'banana', 'cherry']
    ```

    Notes
    -----
    - **Order Preservation**: This function does **not** preserve the original 
      order of elements. If maintaining order is essential, consider using 
      alternative methods such as `dict.fromkeys` or `collections.OrderedDict`.
    - **Hashable Elements**: All elements in the input list must be hashable. 
      Unhashable elements (e.g., lists, dictionaries) will raise a `TypeError`.

    See Also
    --------
    dict.fromkeys : Create a dictionary with keys from the input list, 
        preserving order (Python 3.7+).
    collections.OrderedDict : Ordered dictionary for maintaining element 
        order (pre-Python 3.7).
    """
    return list(set(lst))

def unzip_dict_items(input_dict: Dict[Any, Any]) -> Tuple[List[Any], List[Any]]:
    """
    Unzip a dictionary into two lists containing its keys and values, respectively.

    Parameters:
    ---------------
    input_dict : Dict[Any, Any]
        The dictionary from which to extract keys and values.

    Returns:
    ---------------
    Tuple[List[Any], List[Any]]
        A tuple containing two lists: the first with keys and the second with values from the dictionary.

    Examples:
    ---------------
    >>> test_dict = {'a': 1, 'b': 2, 'c': 3}
    >>> keys, values = unzip_dict_items(test_dict)
    >>> print(keys)
    ['a', 'b', 'c']
    >>> print(values)
    [1, 2, 3]
    """
    keys, values = zip(*input_dict.items()) if input_dict else ([], [])
    return list(keys), list(values)

