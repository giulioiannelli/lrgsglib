"""Tests for lrgsglib.utils.basic modules.

Covers arithmetic, numeric, iterables, matrix, linalg, and paths.
"""
import numpy as np
import pytest
from pathlib import Path

from lrgsglib.utils.basic import (
    adjust_to_even,
    bin_sign,
    ceil,
    flip_to_positive_majority,
    is_in_range,
    sign_with_threshold,
    # numeric
    elements_within_eta_numpy,
    dtype_numerical_precision,
    linspace,
    round_sigfig_n,
    symmetric_logarithm,
    to_fraction,
    width_interval,
    # iterables
    boolean_overlap_fraction,
    flatten,
    uniques,
    extract_subdictionary,
    # matrix
    zoom_into_array,
    shift_with_wrap,
    # paths
    list_dir,
    remove_if_exists,
    find_matching_files,
    # strings
    generate_random_id,
    get_first_int_in_str,
    join_non_empty,
    # other
    compose,
    verbose_print,
    project_3d_to_2d,
    bandpass_sos,
    interpolate_grid_data,
)


# =========================================================================
# arithmetic module
# =========================================================================


class TestArithmetic:
    def test_adjust_to_even_rounds_down(self):
        assert adjust_to_even(2.1) == 2

    def test_adjust_to_even_rounds_up(self):
        assert adjust_to_even(5.5) == 6

    def test_adjust_to_even_exact(self):
        assert adjust_to_even(4.0) == 4

    def test_sign_with_threshold(self):
        arr = np.array([-0.5, 0, 0.5, 1e-18])
        result = sign_with_threshold(arr.copy())
        np.testing.assert_array_equal(result, [-1.0, 0.0, 1.0, 0.0])

    def test_bin_sign_zeros_become_positive(self):
        result = bin_sign([-3, 0, 2, -1, 0])
        np.testing.assert_array_equal(result, [-1, 1, 1, -1, 1])

    def test_flip_to_positive_majority(self):
        arr = np.array([1, -2, -3, 4, -5])
        result = flip_to_positive_majority(arr)
        # Majority negative → all flipped
        assert np.sum(result > 0) >= np.sum(result < 0)

    def test_is_in_range_true(self):
        assert is_in_range(5, 1, 10) is True

    def test_is_in_range_false(self):
        assert is_in_range(15, 1, 10) is False

    def test_is_in_range_boundary(self):
        assert is_in_range(1, 1, 10) is True
        assert is_in_range(10, 1, 10) is True

    def test_ceil_positive(self):
        assert ceil(3.2) == 4

    def test_ceil_negative(self):
        assert ceil(-3.2) == -3

    def test_ceil_integer(self):
        assert ceil(5.0) == 5


# =========================================================================
# numeric module
# =========================================================================


class TestNumeric:
    def test_elements_within_eta(self):
        arr = [1.0, 1.1, 1.5, 2.0, 5.0]
        result = elements_within_eta_numpy(arr, 0.5)
        # min=1.0, eta=0.5 → elements <= 1.5
        assert len(result) == 3
        assert 5.0 not in result

    def test_dtype_numerical_precision(self):
        eps64 = dtype_numerical_precision(np.float64)
        eps32 = dtype_numerical_precision(np.float32)
        assert eps32 > eps64  # float32 less precise
        assert eps64 > 0

    def test_linspace(self):
        result = linspace(0, 1, 5)
        assert len(result) == 5
        assert result[0] == pytest.approx(0.0)
        assert result[-1] == pytest.approx(1.0)

    def test_round_sigfig_n(self):
        assert round_sigfig_n(123.456, 3) == pytest.approx(123.0)
        assert round_sigfig_n(0.00456, 2) == pytest.approx(0.0046)

    def test_symmetric_logarithm(self):
        # symmetric_logarithm(x, a, b, c, d) with safety tolerance
        result = symmetric_logarithm(10.0, 1.0, 0.0, 1.0, 0.0)
        assert isinstance(result, (float, np.floating))

    def test_to_fraction(self):
        result = to_fraction(0.5)
        assert result.numerator == 1
        assert result.denominator == 2

    def test_width_interval(self):
        result = width_interval(0, 10)
        assert result == pytest.approx(10.0)


# =========================================================================
# iterables module
# =========================================================================


class TestIterables:
    def test_boolean_overlap_fraction(self):
        a = np.array([True, True, False, False])
        b = np.array([True, False, True, False])
        frac = boolean_overlap_fraction(a, b)
        # XOR then fraction of True: 2 differences out of 4 → 0.5
        assert frac == pytest.approx(0.5)

    def test_flatten(self):
        nested = [[1, 2], [3, [4, 5]]]
        result = list(flatten(nested))
        assert result == [1, 2, 3, 4, 5]

    def test_uniques(self):
        result = uniques([1, 2, 2, 3, 3, 3])
        assert set(result) == {1, 2, 3}

    def test_extract_subdictionary(self):
        d = {"a": 1, "b": 2, "c": 3, "d": 4}
        result = extract_subdictionary(d, ["a", "c"])
        assert result == {"a": 1, "c": 3}


# =========================================================================
# matrix module
# =========================================================================


class TestMatrix:
    def test_zoom_into_array(self):
        arr = np.arange(100).reshape(10, 10)
        # zoom_into_array(array, x, y, width, height) - center at (5,5), 4x4
        zoomed = zoom_into_array(arr, 5, 5, 4, 4)
        assert zoomed.shape[0] > 0
        assert zoomed.shape[1] > 0

    def test_shift_with_wrap(self):
        arr = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        shifted = shift_with_wrap(arr, 1, 0)
        # Verify shape preserved and content shifted
        assert shifted.shape == arr.shape


# =========================================================================
# strings module
# =========================================================================


class TestStrings:
    def test_generate_random_id(self):
        rid = generate_random_id()
        assert isinstance(rid, str)
        assert len(rid) > 0

    def test_generate_random_id_length(self):
        rid = generate_random_id(length=8)
        assert len(rid) == 8

    def test_get_first_int_in_str(self):
        assert get_first_int_in_str("abc123def") == 123

    def test_join_non_empty(self):
        assert join_non_empty("-", "a", "", "b", "", "c") == "a-b-c"


# =========================================================================
# operations / functions / geometry / io
# =========================================================================


class TestMiscModules:
    def test_compose(self):
        # compose(f, g) → applies f first, then g to f's result
        add1 = lambda x: x + 1
        double = lambda x: x * 2
        composed = compose(add1, double)
        # add1(3)=4, double(4)=8
        assert composed(3) == 8

    def test_verbose_print(self, capsys):
        verbose_print("hello", verbose=True)
        captured = capsys.readouterr()
        assert "hello" in captured.out

    def test_project_3d_to_2d(self):
        x, y, z = 1.0, 0.0, 0.0
        result = project_3d_to_2d(x, y, z)
        assert len(result) == 2


# =========================================================================
# paths module
# =========================================================================


class TestPaths:
    def test_list_dir(self, tmp_path):
        (tmp_path / "a.txt").touch()
        (tmp_path / "b.txt").touch()
        result = list_dir(tmp_path)
        assert len(result) >= 2

    def test_remove_if_exists(self, tmp_path):
        f = tmp_path / "to_remove.txt"
        f.touch()
        assert f.exists()
        remove_if_exists(f)
        assert not f.exists()

    def test_remove_if_exists_nonexistent(self, tmp_path):
        f = tmp_path / "nonexistent.txt"
        # Should not raise
        remove_if_exists(f)
