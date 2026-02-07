"""
Unit tests for SlaplSpect kernel module.

Tests the generic spectral analysis framework used by L2D_SlaplSpect,
MCG_SlaplSpect, and other spectral analysis modules.
"""

import pytest
import numpy as np
import pickle as pk
from pathlib import Path
from collections import Counter
from types import SimpleNamespace

import sys
# sys.path needed for kernels imports (not installed packages)
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from kernels.SlaplSpect import (
    save_data,
    load_existing_data,
    check_existing_and_needed_averages,
    create_initial_bin_counter,
    build_eigval_fname_base,
    build_eigvec_fname_base,
    process_eigen_distribution,
)


class TestSaveLoadData:
    """Test data persistence functions."""

    def test_save_data_counter(self, tmp_path):
        """Test saving Counter data."""
        args = SimpleNamespace(verbose=False)
        data = Counter([1, 2, 2, 3, 3, 3])
        fname = tmp_path / "test_counter.pkl"

        save_data(args, data, fname)

        assert fname.exists()
        with open(fname, "rb") as f:
            loaded = pk.load(f)
        assert loaded == data

    def test_save_data_list(self, tmp_path):
        """Test saving list of arrays."""
        args = SimpleNamespace(verbose=False)
        data = [np.array([1, 2, 3]), np.array([4, 5, 6])]
        fname = tmp_path / "test_list.pkl"

        save_data(args, data, fname)

        assert fname.exists()
        with open(fname, "rb") as f:
            loaded = pk.load(f)
        assert len(loaded) == len(data)
        np.testing.assert_array_equal(loaded[0], data[0])

    def test_load_existing_data_eigval(self, tmp_path):
        """Test loading eigval_dist checkpoint."""
        # Create fake checkpoint
        bin_counter = Counter([10, 20, 20, 30, 30, 30])
        fname = tmp_path / "checkpoint.pkl"
        with open(fname, "wb") as f:
            pk.dump(bin_counter, f)

        # Load it back
        loaded_counter, initial_data = load_existing_data(
            fname, "eigval_dist", verbose=False
        )

        assert loaded_counter == bin_counter
        assert sorted(initial_data) == [10, 20, 20, 30, 30, 30]

    def test_load_existing_data_eigvec(self, tmp_path):
        """Test loading eigvec_dist checkpoint (list of Counters)."""
        # Create fake checkpoint with 3 eigenvectors
        bin_counter = [
            Counter([1, 1, 2]),
            Counter([5, 5, 5]),
            Counter([10])
        ]
        fname = tmp_path / "checkpoint.pkl"
        with open(fname, "wb") as f:
            pk.dump(bin_counter, f)

        # Load it back
        loaded_counter, initial_data = load_existing_data(
            fname, "eigvec_dist", verbose=False
        )

        assert len(loaded_counter) == 3
        assert loaded_counter[0] == bin_counter[0]
        # initial_data should be flattened from all counters
        assert sorted(initial_data) == [1, 1, 2, 5, 5, 5, 10]


class TestCheckpointDetection:
    """Test checkpoint file detection logic."""

    def test_no_checkpoints(self, tmp_path):
        """Test when no checkpoint files exist."""
        navg, path = check_existing_and_needed_averages(
            tmp_path, "test_base", verbose=False
        )

        assert navg == 0
        assert path == tmp_path / "test_base_na=0.pkl"

    def test_single_checkpoint(self, tmp_path):
        """Test finding a single checkpoint."""
        # Create a fake checkpoint file
        fname = tmp_path / "test_base_na=100.pkl"
        fname.touch()

        navg, path = check_existing_and_needed_averages(
            tmp_path, "test_base", verbose=False
        )

        assert navg == 100
        assert path == tmp_path / "test_base_na=100.pkl"

    def test_multiple_checkpoints(self, tmp_path):
        """Test finding the latest checkpoint among multiple files."""
        # Create multiple checkpoint files
        (tmp_path / "test_base_na=10.pkl").touch()
        (tmp_path / "test_base_na=50.pkl").touch()
        (tmp_path / "test_base_na=100.pkl").touch()

        navg, path = check_existing_and_needed_averages(
            tmp_path, "test_base", verbose=False
        )

        assert navg == 100  # Should find the maximum
        assert path == tmp_path / "test_base_na=100.pkl"

    def test_nonexistent_directory(self, tmp_path):
        """Test when directory doesn't exist."""
        fake_dir = tmp_path / "nonexistent"

        navg, path = check_existing_and_needed_averages(
            fake_dir, "test_base", verbose=False
        )

        assert navg == 0


class TestBinCounterCreation:
    """Test bin counter initialization."""

    def test_create_eigval_counter(self):
        """Test creating single Counter for eigval_dist."""
        counter = create_initial_bin_counter("eigval_dist")

        assert isinstance(counter, Counter)
        assert len(counter) == 0

    def test_create_eigvec_counter(self):
        """Test creating list of Counters for eigvec_dist."""
        counter = create_initial_bin_counter("eigvec_dist", howmany=5)

        assert isinstance(counter, list)
        assert len(counter) == 5
        assert all(isinstance(c, Counter) for c in counter)
        assert all(len(c) == 0 for c in counter)


class TestFilenameBuilders:
    """Test filename construction functions."""

    def test_build_eigval_fname_lattice(self):
        """Test eigenvalue filename for lattices."""
        fname = build_eigval_fname_base("", 0.05, "square")
        assert fname == "dist_eigval_p=0.05_square"

    def test_build_eigval_fname_mcg(self):
        """Test eigenvalue filename for MCG."""
        fname = build_eigval_fname_base("mc_", 0.1)
        assert fname == "mc_dist_eigval_p=0.1"

    def test_build_eigvec_fname_lattice(self):
        """Test eigenvector filename for lattices."""
        fname = build_eigvec_fname_base("", 10, 0.05, "smallest")
        assert fname == "dist10_p=0.05_smallest"

    def test_build_eigvec_fname_mcg(self):
        """Test eigenvector filename for MCG."""
        fname = build_eigvec_fname_base("mc_", 5, 0.1, "largest")
        assert fname == "mc_dist5_p=0.1_largest"

    def test_build_fname_precision(self):
        """Test that floating point formatting works correctly."""
        # Test various p values
        assert "p=0.001" in build_eigval_fname_base("", 0.001)
        assert "p=0.5" in build_eigval_fname_base("", 0.5)
        assert "p=1" in build_eigval_fname_base("", 1.0)


class TestProcessEigenDistribution:
    """Integration test for process_eigen_distribution."""

    def test_simple_workflow(self, tmp_path):
        """Test a simple end-to-end workflow."""
        # Mock arguments
        args = SimpleNamespace(
            verbose=False,
            mode="eigval_dist",
            number_of_averages=10,
            save_frequency=5,
            bins_count=20
        )

        # Mock initial data function
        def initial_fn(args):
            return np.array([1.0, 2.0, 3.0, 4.0, 5.0])

        # Mock update function (just returns the counter)
        def update_fn(batch_size, bins, bin_centers, counter, args):
            # Simulate adding some data
            from lrgsglib.config.funcs import bin_eigenvalues
            fake_data = np.random.rand(batch_size * 10)
            counter.update(bin_eigenvalues(fake_data, bins, bin_centers))
            return counter

        # Mock path function
        def get_path_fn(args):
            return tmp_path

        # Run the process
        process_eigen_distribution(
            "test_eigval",
            initial_fn,
            update_fn,
            save_data,
            args,
            get_path_fn
        )

        # Check that final file exists
        final_file = tmp_path / "test_eigval_na=10.pkl"
        assert final_file.exists()

        # Load and verify it's a Counter
        with open(final_file, "rb") as f:
            result = pk.load(f)
        assert isinstance(result, Counter)

    def test_resume_from_checkpoint(self, tmp_path):
        """Test resuming from an existing checkpoint."""
        # Create an existing checkpoint
        existing_counter = Counter([1, 1, 2, 2, 2, 3])
        checkpoint_file = tmp_path / "test_resume_na=5.pkl"
        with open(checkpoint_file, "wb") as f:
            pk.dump(existing_counter, f)

        # Mock arguments requesting 10 total averages
        args = SimpleNamespace(
            verbose=False,
            mode="eigval_dist",
            number_of_averages=10,
            save_frequency=5,
            bins_count=20
        )

        # Mock functions
        def initial_fn(args):
            # This should NOT be called since we have a checkpoint
            raise AssertionError("initial_fn should not be called when resuming")

        def update_fn(batch_size, bins, bin_centers, counter, args):
            # Should only process remaining 5 averages
            assert batch_size == 5
            return counter

        def get_path_fn(args):
            return tmp_path

        # Run the process
        process_eigen_distribution(
            "test_resume",
            initial_fn,
            update_fn,
            save_data,
            args,
            get_path_fn
        )

        # Check that final file exists
        final_file = tmp_path / "test_resume_na=10.pkl"
        assert final_file.exists()

        # Old checkpoint should be removed after completion
        assert not checkpoint_file.exists()


class TestIntegration:
    """Integration tests verifying the refactoring preserves behavior."""

    def test_eigval_dist_workflow(self, tmp_path):
        """Test complete eigval_dist workflow."""
        args = SimpleNamespace(
            verbose=False,
            mode="eigval_dist",
            number_of_averages=20,
            save_frequency=10,
            bins_count=50
        )

        call_count = {"initial": 0, "update": 0}

        def initial_fn(args):
            call_count["initial"] += 1
            return np.random.randn(100)

        def update_fn(batch_size, bins, bin_centers, counter, args):
            call_count["update"] += 1
            from lrgsglib.config.funcs import bin_eigenvalues
            data = np.random.randn(batch_size * 100)
            counter.update(bin_eigenvalues(data, bins, bin_centers))
            return counter

        def get_path_fn(args):
            return tmp_path

        process_eigen_distribution(
            "integration_test",
            initial_fn,
            update_fn,
            save_data,
            args,
            get_path_fn
        )

        # Verify function calls
        assert call_count["initial"] == 1  # Called once for binning
        assert call_count["update"] == 2   # Called twice (2 periods of 10)

        # Verify final file
        final_file = tmp_path / "integration_test_na=20.pkl"
        assert final_file.exists()

    def test_eigvec_dist_workflow(self, tmp_path):
        """Test complete eigvec_dist workflow with list of Counters."""
        args = SimpleNamespace(
            verbose=False,
            mode="eigvec_dist",
            number_of_averages=15,
            save_frequency=5,
            bins_count=30,
            howmany=3
        )

        def initial_fn(args):
            # Return eigenvector components (3 eigenvectors)
            return np.random.rand(3, 100)

        def update_fn(batch_size, bins, bin_centers, counter, args):
            from lrgsglib.config.funcs import bin_eigenvalues
            # Update each eigenvector's counter
            for i in range(args.howmany):
                data = np.random.rand(batch_size * 100)
                counter[i].update(bin_eigenvalues(data, bins, bin_centers))
            return counter

        def get_path_fn(args):
            return tmp_path

        process_eigen_distribution(
            "eigvec_test",
            initial_fn,
            update_fn,
            save_data,
            args,
            get_path_fn
        )

        # Verify final file
        final_file = tmp_path / "eigvec_test_na=15.pkl"
        assert final_file.exists()

        # Load and verify structure
        with open(final_file, "rb") as f:
            result = pk.load(f)

        assert isinstance(result, list)
        assert len(result) == 3
        assert all(isinstance(c, Counter) for c in result)


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])
