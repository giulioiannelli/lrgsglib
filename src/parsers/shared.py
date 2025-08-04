import argparse
import ast
import numpy as np

def parse_arguments(parser: argparse.ArgumentParser) -> argparse.Namespace:
    args = parser.parse_args()
    return args

def parsDict_get(parsDict, var, key):
    return parsDict[var].get(key, None)

def parse_linspace(value):
    """Parses a tuple string and returns a linspace array."""
    try:
        start, stop, num = ast.literal_eval(value)  # Safe parsing of tuples
        return np.linspace(float(start), float(stop), int(num))
    except (ValueError, SyntaxError, TypeError):
        raise argparse.ArgumentTypeError("Invalid linspace format. Use (start, stop, num)")

def parse_multiple_linspace(value):
    """Parses multiple linspace tuples separated by ';' and concatenates them."""
    segments = value.split(';')  # Split multiple tuples
    arrays = [parse_linspace(seg) for seg in segments]  # Convert each to linspace
    return np.concatenate(arrays)  # Merge all into one array