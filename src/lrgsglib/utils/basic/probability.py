from .common import *
#
__all__ = [
    'binder_cumulant',
    'coarsen_bins_with_padding',
    'linear_binning_hist',
    'create_symmetric_log_bins',
    'log_binning',
    'neglog_binning',
    'symlog_binning',
    'marchenko_pastur',
    'update_mean_var',
    'update_mean_m2',
    'update_mean_var_chunk'
]
#
def binder_cumulant(data):
    """
    Calculate the Binder cumulant for a set of data.

    Parameters:
    -----------
    data : array_like
        Input data representing measurements of the order parameter.

    Returns:
    --------
    float
        Binder cumulant of the data.

    Notes:
    ------
    The Binder cumulant is a statistical measure used in the analysis of 
    phase transitions in statistical physics. It is defined as 
    1 - (mean(m^4) / (3 * mean(m^2)^2)), where m is the order parameter.
    The Binder cumulant is particularly useful in finite-size scaling 
    analysis as it is dimensionless and often has a universal value at 
    the critical point for systems within the same universality class.

    Examples:
    ---------
    >>> data = np.random.normal(size=1000)
    >>> print(binder_cumulant(data))

    >>> uniform_data = np.random.uniform(-1, 1, size=1000)
    >>> print(binder_cumulant(uniform_data))
    """
    m2 = np.mean(data**2)
    m4 = np.mean(data**4)
    U4 = 1 - m4 / (3 * m2**2)
    return U4
#
def coarsen_bins_with_padding(x: NDArray, y: NDArray, factor: int) -> Tuple[NDArray, NDArray]:
    """
    Coarsen binned data by a given factor with padding.

    This function groups consecutive bins by the specified factor. The new bin centers are the
    mean of each group, and the new bin values are the sum of each group. If the length of the input
    arrays is not divisible by the factor, they are padded accordingly.

    Parameters
    ----------
    x : NDArray
        1D array of bin centers.
    y : NDArray
        1D array of bin values.
    factor : int
        Factor by which to reduce the number of bins.

    Returns
    -------
    Tuple[NDArray, NDArray]
        Tuple containing:
        - new_x: 1D array of coarsened bin centers.
        - new_y: 1D array of coarsened bin values.

    Notes
    -----
    If len(x) is not a multiple of factor, x is padded using its edge value and y is padded with zeros.
    """
    remainder = len(x) % factor
    if remainder != 0:
        pad_width = factor - remainder
        x = np.pad(x, (0, pad_width), mode='edge')
        y = np.pad(y, (0, pad_width), mode='constant', constant_values=0)
    
    new_x = x.reshape(-1, factor).mean(axis=1)
    new_y = y.reshape(-1, factor).sum(axis=1)
    return new_x, new_y
#
def linear_binning_hist(
    data: NDArray,
    bins_count: int = 20,
    *,
    include_counts: bool = False
) -> Union[Tuple[NDArray, NDArray], Tuple[NDArray, NDArray, NDArray]]:
    """
    Create linearly spaced bins over the exact data range, compute centers—and
    optionally counts—in one function.

    Parameters
    ----------
    data : NDArray
        Input data array.
    bins_count : int, default 20
        Number of equal-width bins.
    include_counts : bool, keyword-only, default True
        If True, return histogram counts as well.

    Returns
    -------
    (bin_edges, bin_centers[, counts])
        - bin_edges : ndarray of length (bins_count + 1)
        - bin_centers : ndarray of length bins_count
        - counts : ndarray of length bins_count (only if `include_counts`)

    """
    bin_edges = np.linspace(data.min(), data.max(), bins_count + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    if include_counts:
        counts, _ = np.histogram(data, bins=bin_edges)
        return bin_edges, bin_centers, counts
    return bin_edges, bin_centers
#

def create_symmetric_log_bins(data: NDArray, num_bins: int, magnitude_increment: int = 2) -> Tuple[NDArray, NDArray]:
    """
    Creates symmetric logarithmic bins based on the specified data, including bin centers.

    Parameters
    ----------
    data : NDArray
        The input data array for which symmetric logarithmic bins are created.
    num_bins : int
        The total number of bins (should be even to allow symmetry).
    magnitude_increment : int, optional
        Increment to control the logarithmic scale extension. Default is 2.

    Returns
    -------
    Tuple[NDArray, NDArray]
        bins : NDArray
            The edges of the bins.
        bin_centers : NDArray
            The center values of each bin.
    """
    min_value = np.min(np.abs(data[data != 0]))
    max_value = np.max(np.abs(data))

    # Creating positive and negative symmetric logarithmic bins
    positive_bins = np.logspace(np.log10(min_value) - magnitude_increment, np.log10(max_value) + magnitude_increment, int(num_bins // 2) + 1)
    negative_bins = -np.flip(positive_bins[:-1])
    
    bins = np.concatenate((negative_bins, [0], positive_bins))
    bin_centers = (bins[:-1] + bins[1:]) / 2

    return bins, bin_centers
#
def log_binning(
    data: NDArray[np.floating],
    bins_count: int = 20
) -> Tuple[NDArray[np.floating], NDArray[np.integer], NDArray[np.floating]]:
    """
    Compute a histogram with logarithmically spaced bins.

    Parameters
    ----------
    data : NDArray[np.floating]
        One-dimensional array of positive values.
    bins_count : int, optional
        Number of log-spaced bins (default is 20).

    Returns
    -------
    bin_centers : NDArray[np.floating]
        Geometric centers of each bin.
    counts : NDArray[int]
        Number of samples falling into each bin.
    widths : NDArray[np.floating]
        Width of each bin (edge_i+1 − edge_i).

    Raises
    ------
    ValueError
        If any entry in `data` is non-positive, since log-binning
        requires strictly positive values.
    """
    if data.ndim != 1:
        raise ValueError("`data` must be a 1D array.")
    if np.any(data <= 0):
        raise ValueError("All `data` values must be > 0 for log binning.")

    log_min = np.log10(data.min())
    log_max = np.log10(data.max())
    # +1 so that we get bins_count intervals, not points
    bin_edges = np.logspace(log_min, log_max, bins_count + 1, base=10)
    counts, _ = np.histogram(data, bins=bin_edges)
    # geometric mean is a better “center” for log spacing
    bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])
    bin_widths = bin_edges[1:] - bin_edges[:-1]

    return bin_centers, counts, bin_widths
#
def neglog_binning(
    data: Union[Sequence[float], NDArray[np.floating]],
    bins_count: int = 20
) -> Tuple[NDArray[np.floating], NDArray[np.integer], NDArray[np.floating]]:
    """
    Compute a histogram with bins equally spaced in the negative log10 scale.

    Parameters
    ----------
    data : array-like of float
        One-dimensional collection of strictly negative values.
    bins_count : int, optional
        Number of log-spaced bins (default: 20).

    Returns
    -------
    bin_centers : ndarray of float
        Geometric centers of bins on the negative axis.
    counts : ndarray of int
        Number of samples in each bin.
    widths : ndarray of float
        Width of each bin.

    Raises
    ------
    ValueError
        If `data` is empty, not one-dimensional, or contains any ≥ 0 values.
    """
    arr = np.asarray(data)
    if arr.ndim != 1:
        raise ValueError("`data` must be one-dimensional.")
    if arr.size == 0:
        raise ValueError("`data` must not be empty.")
    if np.any(arr >= 0):
        raise ValueError("All `data` values must be negative for negative log binning.")

    abs_data = np.abs(arr)
    log_min = np.floor(np.log10(abs_data.min()))
    log_max = np.ceil(np.log10(abs_data.max()))

    # Create positive log-spaced edges, then convert to negative ascending edges
    pos_edges = np.logspace(log_min, log_max, bins_count + 1, base=10)
    bin_edges = -pos_edges[::-1]

    counts, _ = np.histogram(arr, bins=bin_edges)

    # Geometric centers on negative axis
    bin_centers = -np.sqrt(pos_edges[:-1] * pos_edges[1:])[::-1]
    bin_widths = bin_edges[1:] - bin_edges[:-1]

    return bin_centers, counts, bin_widths

def symlog_binning(
    full_data: NDArray[np.floating],
    bins_count: int = 20
) -> Tuple[
    Optional[Tuple[NDArray[np.floating], NDArray[np.integer], NDArray[np.floating]]],
    Optional[Tuple[NDArray[np.floating], NDArray[np.integer], NDArray[np.floating]]]
]:
    """
    Perform symmetric log-histogramming on positive and negative values separately.

    Splits `full_data` into positive (>0) and negative (<0) parts, then calls
    `log_binning` on the positives and `neglog_binning` on the negatives.

    Parameters
    ----------
    full_data : array_like of float
        1D array containing both positive and negative samples.
    bins_count : int, optional
        Number of log-spaced bins for each half (default: 20).

    Returns
    -------
    (pos_result, neg_result) : tuple
        - pos_result : (centers_pos, counts_pos, widths_pos) or None if no positives
        - neg_result : (centers_neg, counts_neg, widths_neg) or None if no negatives

    Raises
    ------
    ValueError
        If `full_data` is not one-dimensional.
    """
    arr = np.asarray(full_data, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"`full_data` must be 1D, got shape {arr.shape}")

    # Positive side
    pos = arr[arr > 0]
    pos_result = log_binning(pos, bins_count) if pos.size else None

    # Negative side
    neg = arr[arr < 0]
    neg_result = neglog_binning(neg, bins_count) if neg.size else None

    return pos_result, neg_result
#
def marchenko_pastur(l, g):
    """
    Compute the Marchenko-Pastur distribution density.

    Parameters
    ----------
    l : array_like
        Input eigenvalue(s) at which the density is evaluated.
    g : float
        Aspect ratio parameter, typically defined as N/M or M/N.

    Returns
    -------
    ndarray
        The Marchenko-Pastur density evaluated at `l`.

    Notes
    -----
    The density is defined as:

        p(l) = (1 / (2 * π * g * l)) * sqrt( max((1+√g)² - l, 0) * max(l - (1-√g)², 0) )

    Reference
    ---------
    Marchenko, V. A. and Pastur, L. A. (1967). "Distribution of eigenvalues for some sets 
    of random matrices." Mathematics of the USSR-Sbornik, 1(4), 457-483.
    """
    def m0(a):
        """Compute the element-wise maximum of the array `a` and 0."""
        return np.maximum(a, 0)

    g_plus = (1 + np.sqrt(g))**2
    g_minus = (1 - np.sqrt(g))**2

    density = np.sqrt(m0(g_plus - l) * m0(l - g_minus)) / (2 * np.pi * g * l)
    return density
#
def update_mean_var(
    mean: NDArray,
    var: NDArray,
    count: int,
    sample: NDArray
) -> Tuple[NDArray, NDArray, int]:
    """
    Online update of mean and variance (Welford's algorithm, vectorized).

    Parameters
    ----------
    mean : NDArray of shape (N,)
        Current mean of each variable over `count` samples.
    var : NDArray of shape (N,)
        Current variance of each variable (population variance) over `count` samples.
    count : int
        Number of samples used to compute `mean` and `var`.
    sample : NDArray of shape (N,)
        New realization to include.

    Returns
    -------
    new_mean : NDArray of shape (N,)
    new_var : NDArray of shape (N,)
    new_count : int
    """
    new_count = count + 1
    delta = sample - mean
    new_mean = mean + delta / new_count
    # update population variance
    new_var = (var * count + delta * (sample - new_mean)) / new_count
    return new_mean, new_var, new_count
#
def update_mean_m2(
    mean: NDArray[np.float64],
    M2:   NDArray[np.float64],
    count: int,
    sample: NDArray[np.float64]
) -> Tuple[NDArray[np.float64], NDArray[np.float64], int]:
    """
    Welford’s online update, vectorized.
    mean_new, M2_new, count_new = update_mean_m2(mean, M2, count, sample)
    """
    count += 1
    delta  = sample - mean
    mean  += delta / count
    delta2 = sample - mean
    M2    += delta * delta2
    return mean, M2, count
#
def update_mean_var_chunk(
    mean: NDArray,
    var: NDArray,
    count: int,
    chunk: NDArray
) -> Tuple[NDArray, NDArray, int]:
    """
    Update mean and variance with a batch of samples.

    Parameters
    ----------
    mean : NDArray, shape (N,)
        Current mean.
    var : NDArray, shape (N,)
        Current population variance.
    count : int
        Number of samples in current statistics.
    chunk : NDArray, shape (M, N)
        New samples batch.

    Returns
    -------
    new_mean : NDArray, shape (N,)
    new_var : NDArray, shape (N,)
    new_count : int
    """
    m2 = chunk.shape[0]
    if m2 == 0:
        return mean, var, count

    # compute chunk statistics
    chunk_mean = chunk.mean(axis=0)
    chunk_var = chunk.var(axis=0)  # population variance

    total_count = count + m2
    delta = chunk_mean - mean

    combined_mean = (count * mean + m2 * chunk_mean) / total_count

    # combined second moment
    ss1 = var * count
    ss2 = chunk_var * m2
    ss_delta = count * m2 * delta**2 / total_count

    combined_var = (ss1 + ss2 + ss_delta) / total_count

    return combined_mean, combined_var, total_count