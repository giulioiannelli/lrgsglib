"""Utilities for reading and writing protein structures."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Optional, Tuple, List

import numpy as np
import requests

# Bio.PDB imports for PDB processing
try:
    from Bio.PDB import PDBParser
    BIO_PDB_AVAILABLE = True
except ImportError:  # pragma: no cover - optional dependency
    BIO_PDB_AVAILABLE = False

__all__ = [
    "download_pdb",
    "extract_ca_coordinates",
    "extract_atoms_coordinates",
    "coords_to_pdb_string_with_structure",
]


def download_pdb(pdb_id: str, storage_dir: Optional[Path] = None) -> str:
    """Download PDB file from RCSB with simple caching."""
    if storage_dir is None:
        storage_dir = Path.cwd() / "pdb_cache"
    storage_dir.mkdir(exist_ok=True)

    pdb_file_path = storage_dir / f"{pdb_id.upper()}.pdb"
    if pdb_file_path.exists():
        print(f"  PDB {pdb_id} found in cache: {pdb_file_path}")
        return pdb_file_path.read_text()

    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        pdb_file_path.write_text(response.text)
        print(f"  Downloaded {pdb_id} → {pdb_file_path}")
        return response.text
    raise Exception(f"Failed to download PDB {pdb_id}: HTTP {response.status_code}")


def extract_ca_coordinates(pdb_content: str) -> np.ndarray:
    """Extract CA (alpha carbon) coordinates from a PDB string."""
    if not BIO_PDB_AVAILABLE:
        raise ImportError("Bio.PDB required for PDB parsing. Install with: pip install biopython")

    parser = PDBParser(QUIET=True)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as tmp:
        tmp.write(pdb_content)
        tmp_path = tmp.name
    try:
        structure = parser.get_structure("protein", tmp_path)
        coords = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if "CA" in residue:
                        coords.append(residue["CA"].get_coord())
        return np.array(coords) if coords else np.empty((0, 3))
    finally:
        os.unlink(tmp_path)


def extract_atoms_coordinates(pdb_content: str) -> Tuple[np.ndarray, List[str], List[str], List[int]]:
    """Extract all atom coordinates from PDB content with metadata."""
    if not BIO_PDB_AVAILABLE:
        raise ImportError("Bio.PDB required for PDB parsing. Install with: pip install biopython")

    parser = PDBParser(QUIET=True)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as tmp:
        tmp.write(pdb_content)
        tmp_path = tmp.name
    try:
        structure = parser.get_structure("protein", tmp_path)
        coordinates: List[np.ndarray] = []
        atom_names: List[str] = []
        residue_names: List[str] = []
        residue_numbers: List[int] = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    residue_id = residue.get_id()[1]
                    residue_name = residue.get_resname()
                    for atom in residue:
                        coordinates.append(atom.get_coord())
                        atom_names.append(atom.get_name())
                        residue_names.append(residue_name)
                        residue_numbers.append(residue_id)
        coords = np.array(coordinates) if coordinates else np.empty((0, 3))
        return coords, atom_names, residue_names, residue_numbers
    finally:
        os.unlink(tmp_path)


# ---------------------------------------------------------------------------
# PDB formatting utilities
# ---------------------------------------------------------------------------

def coords_to_pdb_string_with_structure(
    coords: np.ndarray,
    residue_types: Optional[List[str]] = None,
    ss_assignments: Optional[List[str]] = None,
) -> str:
    """Return a simple PDB string from coordinates and secondary structure."""
    if residue_types is None:
        residue_types = ["ALA"] * len(coords)

    if ss_assignments is None:
        ss_assignments = ["H" if i % 8 < 4 else "E" if i % 8 < 6 else "C" for i in range(len(coords))]

    pdb_lines: List[str] = []
    helix_count = 1
    sheet_count = 1
    in_helix = False
    in_sheet = False
    helix_start = None
    sheet_start = None
    for i, ss in enumerate(ss_assignments):
        if ss == "H" and not in_helix:
            helix_start = i + 1
            in_helix = True
        elif ss != "H" and in_helix:
            pdb_lines.append(
                f"HELIX  {helix_count:3d} {helix_count:3d} ALA A{helix_start:4d}  ALA A{i:4d}  1{i-helix_start+1:5d}"
            )
            helix_count += 1
            in_helix = False
        if ss == "E" and not in_sheet:
            sheet_start = i + 1
            in_sheet = True
        elif ss != "E" and in_sheet:
            pdb_lines.append(
                f"SHEET  {sheet_count:3d} A{sheet_count:2d} ALA A{sheet_start:4d}  ALA A{i:4d}  0"
            )
            sheet_count += 1
            in_sheet = False
    if in_helix:
        pdb_lines.append(
            f"HELIX  {helix_count:3d} {helix_count:3d} ALA A{helix_start:4d}  ALA A{len(coords):4d}  1{len(coords)-helix_start+1:5d}"
        )
    if in_sheet:
        pdb_lines.append(
            f"SHEET  {sheet_count:3d} A{sheet_count:2d} ALA A{sheet_start:4d}  ALA A{len(coords):4d}  0"
        )

    CA_C_BOND = 1.525
    C_N_BOND = 1.329
    N_CA_BOND = 1.458
    C_O_BOND = 1.231
    for i, (coord, restype) in enumerate(zip(coords, residue_types)):
        ca_coord = coord
        if i > 0:
            prev_ca = coords[i - 1]
            direction = ca_coord - prev_ca
            direction = direction / np.linalg.norm(direction) if np.linalg.norm(direction) > 0 else np.array([1, 0, 0])
            n_coord = ca_coord - direction * N_CA_BOND
        else:
            n_coord = ca_coord + np.array([-N_CA_BOND, 0, 0])
        if i < len(coords) - 1:
            next_ca = coords[i + 1]
            direction = next_ca - ca_coord
            direction = direction / np.linalg.norm(direction) if np.linalg.norm(direction) > 0 else np.array([1, 0, 0])
            c_coord = ca_coord + direction * CA_C_BOND
        else:
            c_coord = ca_coord + np.array([CA_C_BOND, 0, 0])
        o_offset = np.array([0, C_O_BOND, 0])
        o_coord = c_coord + o_offset
        atom_id = i * 4 + 1
        line_n = (
            f"ATOM  {atom_id:5d}  N   {restype} A{i+1:4d}    "
            f"{n_coord[0]:8.3f}{n_coord[1]:8.3f}{n_coord[2]:8.3f}"
            "  1.00 20.00           N  "
        )
        pdb_lines.append(line_n)
        line_ca = (
            f"ATOM  {atom_id+1:5d}  CA  {restype} A{i+1:4d}    "
            f"{ca_coord[0]:8.3f}{ca_coord[1]:8.3f}{ca_coord[2]:8.3f}"
            "  1.00 20.00           C  "
        )
        pdb_lines.append(line_ca)
        line_c = (
            f"ATOM  {atom_id+2:5d}  C   {restype} A{i+1:4d}    "
            f"{c_coord[0]:8.3f}{c_coord[1]:8.3f}{c_coord[2]:8.3f}"
            "  1.00 20.00           C  "
        )
        pdb_lines.append(line_c)
        line_o = (
            f"ATOM  {atom_id+3:5d}  O   {restype} A{i+1:4d}    "
            f"{o_coord[0]:8.3f}{o_coord[1]:8.3f}{o_coord[2]:8.3f}"
            "  1.00 20.00           O  "
        )
        pdb_lines.append(line_o)
    pdb_lines.append("END")
    return "\n".join(pdb_lines)
