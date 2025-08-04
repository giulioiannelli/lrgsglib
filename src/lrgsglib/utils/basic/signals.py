import warnings

from scipy.signal import butter, sosfiltfilt
from numpy.typing import NDArray

__all__ = [
    "bandpass_sos"
]

def bandpass_sos(
    data: NDArray,
    lowcut: float,
    highcut: float,
    fs: float,
    order: int
) -> NDArray:
    """
    Apply a zero-phase Butterworth bandpass filter using second-order sections (SOS).

    Parameters
    ----------
    data : ndarray
        Input signal of shape (n_channels, n_samples).
    lowcut : float
        Lower cutoff frequency (Hz).
    highcut : float
        Upper cutoff frequency (Hz).
    fs : float
        Sampling frequency of `data` (Hz).
    order : int
        Filter order.

    Returns
    -------
    ndarray
        Bandpass‐filtered signal, same shape as `data`.
        
    Raises
    ------
    ValueError
        If filter parameters are invalid or would cause numerical instabilities.
    """
    # Parameter validation
    nyquist = 0.5 * fs
    
    if lowcut <= 0:
        raise ValueError(f"lowcut must be positive, got {lowcut}")
    if highcut >= nyquist:
        raise ValueError(f"highcut ({highcut}) must be less than Nyquist frequency ({nyquist})")
    if lowcut >= highcut:
        raise ValueError(f"lowcut ({lowcut}) must be less than highcut ({highcut})")
    if order <= 0:
        raise ValueError(f"order must be positive, got {order}")
    
    # Check for potentially problematic filter parameters
    bandwidth_ratio = (highcut - lowcut) / nyquist
    if bandwidth_ratio < 0.01:  # Very narrow band
        warnings.warn(
            f"Very narrow filter bandwidth ({bandwidth_ratio:.4f} of Nyquist). "
            "Consider reducing filter order or increasing bandwidth.",
            RuntimeWarning
        )
    
    if order > 10:  # High order filter
        warnings.warn(
            f"High filter order ({order}) may cause numerical instabilities. "
            "Consider using a lower order filter.",
            RuntimeWarning
        )
    
    # Normalize frequencies with small epsilon to avoid edge cases
    eps = 1e-10
    low_norm = max(lowcut / nyquist, eps)
    high_norm = min(highcut / nyquist, 1.0 - eps)
    
    try:
        sos = butter(
            N=order,
            Wn=[low_norm, high_norm],
            btype='band',
            output='sos'
        )
        return sosfiltfilt(sos, data, axis=1)
    except Exception as e:
        raise ValueError(
            f"Filter design failed with lowcut={lowcut}, highcut={highcut}, "
            f"fs={fs}, order={order}. Error: {str(e)}"
        ) from e