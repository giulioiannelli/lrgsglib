import numpy as np
from sklearn.datasets import fetch_openml


def init_mnist_patterns(digit=None, n_samples=1, threshold=127):
    """Return MNIST patterns converted to bipolar arrays."""
    mnist = fetch_openml("mnist_784", version=1)
    X = mnist.data.astype(np.float32)
    y = mnist.target.astype(np.int64)

    if digit is not None:
        indices = np.where(y == digit)[0]
        if len(indices) == 0:
            raise ValueError(f"No MNIST images found for digit {digit}.")
        X_digit = X.iloc[indices][:n_samples]
    else:
        X_digit = X.iloc[:n_samples]

    patterns = []
    for idx in range(len(X_digit)):
        pattern = X_digit.iloc[idx].values
        pattern = np.where(pattern > threshold, 1, -1)
        patterns.append(pattern)
    return patterns
