from ..common import *
from ..funcs import *
from ..FullyConnected import SignedGraph
from .patterns import init_mnist_patterns


class HofieldNN(SignedGraph):
    """
    HofieldNN is a Hopfield-like neural network that inherits from FullyConnected.
    It uses the fully connected graph structure from FullyConnected and implements
    Hopfield network dynamics.
    """

    def __init__(
        self,
        with_patterns: str = "uniform",
        digit=None,
        n_samples=1,
        threshold=127,
        only_const_mode: bool = False,
        **kwargs,
    ):
        self.only_const_mode = only_const_mode
        # Initialize the underlying fully connected graph
        super().__init__(only_const_mode=only_const_mode, **kwargs)
        if not only_const_mode:
            if with_patterns == "uniform":
                self.patterns = np.random.choice([-1, 1], size=(n_samples, self.N))
            elif with_patterns == "mnist":
                self.patterns = init_mnist_patterns(digit, n_samples, threshold)
            # Calculate the weight matrix using the Hebbian learning rule
            self.compute_hebbian_weights(self.patterns)
            fast_set_weights_from_matrix(self.G, self.weights)

    def __repr__(self):
        return f"HofieldNN(N={self.N}, patterns_stored={len(self.patterns)})"

    def __str__(self):
        return self.__repr__()

    def compute_hebbian_weights(self, patterns):
        """
        Compute the weight matrix using the Hebbian rule.

        Args:
            patterns (list of np.ndarray): The patterns to store.

        Returns:
            np.ndarray: Weight matrix of shape (N, N).
        """
        self.weights = np.zeros((self.N, self.N))
        for pattern in patterns:
            pattern = np.array(pattern).flatten()  # Ensure the pattern is 1D
            self.weights += np.outer(pattern, pattern)
        # Remove self-connections by zeroing the diagonal
        np.fill_diagonal(self.weights, 0)
        self.weights /= len(patterns)

    # def update_state(self, state=None):
    #     """
    #     Update the network state using the Hopfield update rule.

    #     Args:
    #         state (np.ndarray): Optional. The current state vector. If None, uses self.state.

    #     Returns:
    #         np.ndarray: The updated state vector.
    #     """
    #     if state is None:
    #         state = self.state.copy()
    #     new_state = state.copy()
    #     if self.synchronous:
    #         # Synchronous update: update all neurons at once
    #         for i in range(self.N):
    #             net_input = np.dot(self.weights[i], state)
    #             new_state[i] = 1 if net_input >= 0 else -1
    #     else:
    #         # Asynchronous update: update neurons in a random order
    #         indices = np.random.permutation(self.N)
    #         for i in indices:
    #             net_input = np.dot(self.weights[i], new_state)
    #             new_state[i] = 1 if net_input >= 0 else -1
    #     self.state = new_state
    #     return new_state

    # def run_until_convergence(self, initial_state=None, max_steps=100):
    #     """
    #     Run the network until the state converges or until a maximum number of steps is reached.

    #     Args:
    #         initial_state (np.ndarray): Optional initial state. If None, uses self.state.
    #         max_steps (int): Maximum number of update steps.

    #     Returns:
    #         np.ndarray: The converged state vector.
    #     """
    #     if initial_state is not None:
    #         self.state = np.array(initial_state)
    #     for _ in range(max_steps):
    #         previous_state = self.state.copy()
    #         self.update_state()
    #         if np.array_equal(self.state, previous_state):
    #             break
    #     return self.state

    # def energy(self, state=None):
    #     """
    #     Calculate the energy of a given state for the Hopfield network.

    #     Args:
    #         state (np.ndarray): Optional state vector. If None, uses self.state.

    #     Returns:
    #         float: The energy of the network state.
    #     """
    #     if state is None:
    #         state = self.state
    #     return -0.5 * np.dot(state, np.dot(self.weights, state))

    # def init_random_state(self):
    #     """
    #     Initialize the network state randomly with +1 or -1 values for each neuron.

    #     Returns:
    #         np.ndarray: The randomly initialized state.
    #     """
    #     self.state = np.random.choice([1, -1], size=self.N)
    #     return self.state

    # def init_mnist_state(self, digit=None, n_samples=1, threshold=127):
    #     """
    #     Initialize the network state using MNIST digit data.
    #     This method fetches the MNIST dataset, filters for a specific digit if provided,
    #     and converts the selected image(s) to a bipolar pattern (+1/-1) for initialization.
    #     Note: The network size (N) must match the number of pixels in the MNIST images (784).

    #     Args:
    #         digit (int, optional): If provided, only images with this digit are used.
    #         n_samples (int): The number of samples to average over (default is 1).
    #         threshold (int): Pixel threshold for binarization. Pixels greater than this value become 1, else -1.

    #     Returns:
    #         np.ndarray: The initialized state from the MNIST data.
    #     """
    #     # Fetch MNIST dataset from openml
    #     mnist = fetch_openml('mnist_784', version=1)
    #     X = mnist.data.astype(np.float32)
    #     y = mnist.target.astype(np.int64)
    #     if digit is not None:
    #         indices = np.where(y == digit)[0]
    #         if len(indices) == 0:
    #             raise ValueError(f"No MNIST images found for digit {digit}.")
    #         X_digit = X[indices][:n_samples]
    #     else:
    #         X_digit = X[:n_samples]
    #     # If multiple samples are selected, average them; otherwise, use the single image.
    #     pattern = np.mean(X_digit, axis=0) if n_samples > 1 else X_digit[0]
    #     # Convert to bipolar: if pixel > threshold then 1, else -1.
    #     pattern = np.where(pattern > threshold, 1, -1)
    #     if pattern.shape[0] != self.N:
    #         raise ValueError(f"MNIST pattern length {pattern.shape[0]} does not match network size {self.N}")
    #     self.state = pattern
    #     return pattern
