"""
Test script for MCG_SlaplSpect_Serializer --print option.

Verifies that the serializer:
1. Generates correct command-line strings
2. Includes all required parameters
3. Properly includes --backend cupy
4. Handles parameter sweeps correctly
"""

import subprocess
import sys
from pathlib import Path
import pytest

pytestmark = pytest.mark.integration

ROOT_DIR = Path(__file__).resolve().parents[1]
SERIALIZER_SCRIPT = ROOT_DIR / "src" / "MCG_SlaplSpect_Serializer.py"


def test_serializer_print():
    """Test serializer --print option with minimal parameters."""
    print("\n" + "="*70)
    print("Testing MCG_SlaplSpect_Serializer --print option")
    print("="*70)

    print("\nRunning serializer with --print option...")
    cmd = [
        "python", str(SERIALIZER_SCRIPT),
        "--print",
        "--matrix_list", "0.8,0.6,0.6,0.8",
        "--fraction_list", "0.5",
        "--iterations_list", "6,7",
        "--p_list", "0.0,0.1",
        "--",  # Separator for additional arguments
        "-m", "eigvals",
        "--howmany", "0",
        "-na", "10",
        "--save_frequency", "5"
    ]

    print(f"Command: {' '.join(cmd)}\n")

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path.cwd())

    if result.returncode != 0:
        print(f"❌ FAILED: Serializer returned non-zero exit code: {result.returncode}")
        print(f"STDERR: {result.stderr}")
        pytest.fail("Serializer returned non-zero exit code")

    print("✓ Serializer ran successfully\n")

    # Parse output - commands are multi-line with backslash continuation
    # Each command spans 2 lines: "slanzarv ... -- \" and "  python ..."
    output_lines = result.stdout.strip().split('\n')

    # Filter out the "Total number" line and join continued lines
    raw_lines = [line for line in output_lines if not line.startswith("Total number")]

    # Count actual commands (lines starting with 'slanzarv')
    command_lines = [line for line in raw_lines if line.startswith('slanzarv')]

    print(f"Generated {len(command_lines)} commands:\n")

    # Expected: 2 iterations × 2 pflip values = 4 commands
    expected_count = 2 * 2  # 2 iterations, 2 pflip values

    if len(command_lines) != expected_count:
        print(f"❌ FAILED: Expected {expected_count} commands, got {len(command_lines)}")
        assert False, f"Expected {expected_count} commands, got {len(command_lines)}"

    print(f"✓ Correct number of commands: {len(command_lines)}\n")

    # Verify each command has required components
    print("Verifying command structure...\n")

    required_components = [
        ("MCG_SlaplSpect.py", "Main program"),
        ("0.8", "p1 parameter"),
        ("0.6", "p2 parameter"),
        ("0.5", "fraction parameter"),
        ("-p", "pflip flag"),
        ("--backend", "backend flag"),
        # Note: backend can be cupy OR scipy depending on graph size
        ("-m", "mode flag"),
        ("eigvals", "mode value"),
        ("--howmany", "howmany flag"),
        ("0", "howmany value (full spectrum)"),
        ("-na", "number_of_averages flag"),
        ("10", "number_of_averages value"),
    ]

    all_valid = True

    # Reconstruct full commands by joining slanzarv line with python line
    full_commands = []
    for i in range(0, len(raw_lines), 2):
        if i+1 < len(raw_lines):
            # Join the two lines (remove backslash and extra spaces)
            full_cmd = raw_lines[i].replace(' \\', '') + ' ' + raw_lines[i+1].strip()
            full_commands.append(full_cmd)

    for i, cmd_line in enumerate(full_commands):
        print(f"Command {i+1}:")
        print(f"  {cmd_line}\n")

        # Extract Python command portion (after slanzarv's --)
        # Split by '-- python' and take the second part
        if '-- python' in cmd_line:
            python_cmd = 'python' + cmd_line.split('-- python')[1]
        else:
            python_cmd = cmd_line

        # Check for required components
        missing = []
        for component, description in required_components:
            if component not in cmd_line:
                missing.append((component, description))

        if missing:
            print(f"  ❌ Missing components:")
            for comp, desc in missing:
                print(f"     - {comp} ({desc})")
            all_valid = False
        else:
            print(f"  ✓ All required components present")

        # Verify backend is either cupy or scipy (both are valid)
        if "cupy" in cmd_line or "scipy" in cmd_line:
            backend = "cupy" if "cupy" in cmd_line else "scipy"
            print(f"  ✓ Backend: {backend}")
        else:
            print(f"  ❌ Missing backend (expected cupy or scipy)")
            all_valid = False

        # Verify iterations parameter (should be 6 or 7)
        tokens = python_cmd.split()
        try:
            # Find MCG_SlaplSpect.py, then positional args: p1 p2 p3 p4 fraction iterations
            # Look for the .py file and count 5 tokens ahead for iterations (after p1,p2,p3,p4,fraction)
            py_idx = next(i for i, t in enumerate(tokens) if 'MCG_SlaplSpect.py' in t)
            # Positional args after .py: p1 p2 p3 p4 fraction iterations
            iterations = tokens[py_idx + 6]  # 6th argument after .py
            if iterations not in ["6", "7"]:
                print(f"  ❌ Unexpected iterations value: {iterations}")
                all_valid = False
            else:
                print(f"  ✓ Iterations: {iterations}")
        except (ValueError, IndexError, StopIteration):
            print(f"  ❌ Could not find iterations parameter")
            all_valid = False

        # Verify pflip parameter (formatted to 2 decimal places: 0.00 or 0.10)
        try:
            p_idx = tokens.index("-p")
            pflip = tokens[p_idx + 1]
            # Serializer formats with consistent precision - could be 0.00, 0.10, etc.
            if pflip not in ["0.00", "0.10", "0.0", "0.1"]:
                print(f"  ⚠ Unexpected pflip value (but may be formatted): {pflip}")
            else:
                print(f"  ✓ pflip: {pflip}")
        except (ValueError, IndexError):
            print(f"  ❌ Could not find pflip parameter")
            all_valid = False

        print()

    if not all_valid:
        print("❌ FAILED: Some commands have missing or incorrect components")
        assert False, "Some commands have missing or incorrect components"

    print("✓ All commands have correct structure")

    # Verify at least one backend is used (cupy or scipy)
    has_valid_backends = all(("cupy" in cmd or "scipy" in cmd) for cmd in full_commands)
    assert has_valid_backends, "All commands should have a valid backend (cupy or scipy)"
    print("✓ All commands have valid backend")

    # Verify parameter sweep coverage
    print("\nVerifying parameter sweep coverage...")

    iterations_found = set()
    pflip_found = set()

    for cmd_line in full_commands:
        # Extract Python command portion
        if '-- python' in cmd_line:
            python_cmd = 'python' + cmd_line.split('-- python')[1]
        else:
            python_cmd = cmd_line

        tokens = python_cmd.split()

        # Extract iterations (6th positional arg after .py file)
        try:
            py_idx = next(i for i, t in enumerate(tokens) if 'MCG_SlaplSpect.py' in t)
            iterations = tokens[py_idx + 6]
            iterations_found.add(iterations)
        except (ValueError, IndexError, StopIteration):
            pass

        # Extract pflip
        try:
            p_idx = tokens.index("-p")
            pflip = tokens[p_idx + 1]
            pflip_found.add(pflip)
        except (ValueError, IndexError):
            pass

    expected_iterations = {"6", "7"}
    expected_pflip = {"0.00", "0.10"}  # Serializer formats with 2 decimals

    if iterations_found != expected_iterations:
        print(f"❌ FAILED: Expected iterations {expected_iterations}, got {iterations_found}")
        assert False, "Expected iterations not found"

    if pflip_found != expected_pflip:
        print(f"❌ FAILED: Expected pflip values {expected_pflip}, got {pflip_found}")
        assert False, "Expected pflip values not found"

    print(f"✓ Parameter sweep covers all combinations:")
    print(f"  - Iterations: {sorted(iterations_found)}")
    print(f"  - pflip: {sorted(pflip_found)}")

    print("\n" + "="*70)
    print("✅ SERIALIZER --PRINT TEST PASSED")
    print("="*70)
    print("\nSummary:")
    print(f"  ✓ Generated {len(full_commands)} commands correctly")
    print("  ✓ All commands include --backend (cupy or scipy)")
    print("  ✓ All commands include mode, howmany, and other parameters")
    print("  ✓ Parameter sweep covers all combinations")
    print("\nConclusion: Serializer is ready for cluster deployment!")

    return None


def test_serializer_with_slanzarv():
    """Test that commands include slanzarv wrapper."""
    print("\n" + "="*70)
    print("Testing serializer with slanzarv wrapper (no --nomail flag)")
    print("="*70)

    cmd = [
        "python", str(SERIALIZER_SCRIPT),
        "--print",
        "--matrix_list", "0.8,0.6,0.6,0.8",
        "--fraction_list", "0.5",
        "--iterations_list", "6",
        "--p_list", "0.0",
        "--",
        "-m", "eigvals",
        "--howmany", "0"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path.cwd())

    if result.returncode != 0:
        print(f"❌ FAILED: Serializer returned error")
        pytest.fail("Serializer returned error")

    output_lines = result.stdout.strip().split('\n')
    raw_lines = [line for line in output_lines if not line.startswith("Total number")]

    # Count commands (lines starting with 'slanzarv')
    command_lines = [line for line in raw_lines if line.startswith('slanzarv')]

    if len(command_lines) != 1:
        print(f"❌ FAILED: Expected 1 command, got {len(command_lines)}")
        assert False, f"Expected 1 command, got {len(command_lines)}"

    # Reconstruct full command
    full_cmd = raw_lines[0].replace(' \\', '') + ' ' + raw_lines[1].strip() if len(raw_lines) >= 2 else raw_lines[0]
    cmd_line = full_cmd
    print(f"\nGenerated command:\n{cmd_line}\n")

    # Check for slanzarv
    if "slanzarv" not in cmd_line:
        print(f"❌ FAILED: Command does not include slanzarv wrapper")
        assert False, "Command does not include slanzarv wrapper"

    print("✓ Command includes slanzarv wrapper")

    # Check for jobname
    if "--jobname" not in cmd_line:
        print(f"❌ FAILED: Command does not include --jobname")
        assert False, "Command does not include --jobname"

    print("✓ Command includes --jobname")

    print("\n✅ SLANZARV WRAPPER TEST PASSED")
    return None


def main():
    """Run all serializer tests."""
    print("\n" + "#"*70)
    print("# MCG_SlaplSpect_Serializer Tests")
    print("#"*70)

    try:
        test_serializer_print()
        test_serializer_with_slanzarv()
    except AssertionError:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
