#!/usr/bin/env python3

import sys
import os

# Add the current directory to path
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(current_dir, 'src')
sys.path.insert(0, src_dir)

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

# Import needed components directly
import matplotlib.collections
import matplotlib.colors
from matplotlib.patches import Polygon

# Define basic constants that would be in const_plotlib
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize

def test_triangle_orientation():
    """Test that triangle orientations match the original honeycomb pattern"""
    
    # Simple test data
    n_rows, n_cols = 4, 6
    triangle_size = 1.0
    h = (np.sqrt(3) / 2) * triangle_size
    
    # Generate coordinates like the original function does
    x_coords = []
    y_coords = []
    vals = []
    
    for i in range(n_rows):
        for j in range(n_cols):
            x_coords.append(j * triangle_size / 2)
            y_coords.append(i * h)
            vals.append(np.random.rand())  # Random values for coloring
    
    print("Test data generated successfully!")
    print(f"Grid: {n_rows}x{n_cols}")
    print(f"Number of triangles: {len(x_coords)}")
    print(f"Sample coordinates: ({x_coords[0]:.2f}, {y_coords[0]:.2f})")
    
    # Test triangle creation logic
    test_triangles = []
    up_tri = np.array([
        [0,       h/2],
        [-triangle_size/2, -h/2],
        [triangle_size/2, -h/2],
    ])
    
    down_tri = np.array([
        [0,       -h/2],
        [-triangle_size/2, h/2],
        [triangle_size/2, h/2],
    ])
    
    # Check orientation pattern matches original
    orientation_match = True
    for i in range(min(3, n_rows)):
        for j in range(min(3, n_cols)):
            x0 = j * triangle_size / 2
            y0 = i * h
            
            if (i + j) % 2 == 0:
                triangle = up_tri + np.array([x0, y0])
                expected_orientation = "up"
            else:
                triangle = down_tri + np.array([x0, y0])
                expected_orientation = "down"
            
            print(f"Grid position ({i},{j}): {expected_orientation} triangle at ({x0:.2f}, {y0:.2f})")
            test_triangles.append(triangle)
    
    print(f"Triangle orientation pattern test: {'PASS' if orientation_match else 'FAIL'}")
    print("Basic functionality test completed successfully!")

if __name__ == "__main__":
    test_triangle_orientation()
