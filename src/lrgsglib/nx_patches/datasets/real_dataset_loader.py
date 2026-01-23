"""
RealDatasetLoader: Load real-world networks from various sources.

Supports loading from KONECT, SNAP, Network Repository, and local files.
Networks are converted to SignedGraph format with optional sign assignment.

Example
-------
>>> loader = RealDatasetLoader()
>>> G = loader.load("karate")  # Built-in Zachary's karate club
>>> G = loader.load_from_file("my_network.edgelist")
>>> G = loader.load_konect("dolphins")
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union, Literal
import hashlib
import gzip
import io
import os
import tempfile
import warnings

import networkx as nx
import numpy as np

from ..SignedGraph.SignedGraph import SignedGraph


@dataclass
class DatasetInfo:
    """
    Metadata for a real-world network dataset.

    Attributes
    ----------
    name : str
        Dataset name/identifier.
    source : str
        Data source (konect, snap, networkrepository, local, builtin).
    url : str, optional
        Download URL.
    n_nodes : int, optional
        Number of nodes (if known).
    n_edges : int, optional
        Number of edges (if known).
    directed : bool
        Whether the network is directed.
    weighted : bool
        Whether the network has edge weights.
    signed : bool
        Whether the network has signed edges.
    bipartite : bool
        Whether the network is bipartite.
    temporal : bool
        Whether the network has temporal information.
    description : str
        Short description of the dataset.
    citation : str, optional
        Citation for the dataset.
    """

    name: str
    source: str
    url: Optional[str] = None
    n_nodes: Optional[int] = None
    n_edges: Optional[int] = None
    directed: bool = False
    weighted: bool = False
    signed: bool = False
    bipartite: bool = False
    temporal: bool = False
    description: str = ""
    citation: Optional[str] = None
    tags: list[str] = field(default_factory=list)


# Built-in datasets that ship with NetworkX or are easily constructible
BUILTIN_DATASETS: dict[str, DatasetInfo] = {
    "karate": DatasetInfo(
        name="karate",
        source="builtin",
        n_nodes=34,
        n_edges=78,
        directed=False,
        description="Zachary's karate club social network",
        citation="Zachary (1977)",
        tags=["social", "small"],
    ),
    "florentine_families": DatasetInfo(
        name="florentine_families",
        source="builtin",
        n_nodes=15,
        n_edges=20,
        directed=False,
        description="Florentine families marriage network",
        tags=["social", "historical", "small"],
    ),
    "les_miserables": DatasetInfo(
        name="les_miserables",
        source="builtin",
        n_nodes=77,
        n_edges=254,
        directed=False,
        weighted=True,
        description="Character co-appearance in Les Miserables",
        tags=["literary", "small"],
    ),
    "davis_southern_women": DatasetInfo(
        name="davis_southern_women",
        source="builtin",
        n_nodes=32,
        n_edges=89,
        directed=False,
        bipartite=True,
        description="Davis Southern women social events",
        tags=["social", "bipartite", "small"],
    ),
    "football": DatasetInfo(
        name="football",
        source="builtin",
        n_nodes=115,
        n_edges=613,
        directed=False,
        description="American college football network",
        tags=["sports", "small"],
    ),
}


class RealDatasetLoader:
    """
    Load real-world network datasets from various sources.

    Supports loading networks from:
    - Built-in NetworkX datasets
    - Local files (edge list, adjacency list, GraphML, GML, etc.)
    - KONECT (Koblenz Network Collection) - requires internet
    - SNAP (Stanford Network Analysis Project) - requires internet
    - Network Repository - requires internet

    Parameters
    ----------
    cache_dir : str or Path, optional
        Directory for caching downloaded datasets.
        Default: ~/.lrgsglib/datasets
    auto_download : bool, default True
        Whether to automatically download datasets not in cache.
    default_sign : int, default 1
        Default sign to assign to edges without sign information.
    pflip : float, default 0.0
        Probability of flipping edge signs after loading.
    seed : int, optional
        Random seed for reproducible sign flipping.

    Example
    -------
    >>> loader = RealDatasetLoader(cache_dir="/tmp/datasets")
    >>> # Load built-in dataset
    >>> karate = loader.load("karate")
    >>> # Load from local file
    >>> my_net = loader.load_from_file("network.edgelist")
    >>> # Load from KONECT (downloads if needed)
    >>> dolphins = loader.load_konect("dolphins")
    """

    def __init__(
        self,
        cache_dir: Optional[Union[str, Path]] = None,
        auto_download: bool = True,
        default_sign: int = 1,
        pflip: float = 0.0,
        seed: Optional[int] = None,
    ) -> None:
        if cache_dir is None:
            cache_dir = Path.home() / ".lrgsglib" / "datasets"
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        self.auto_download = auto_download
        self.default_sign = default_sign
        self.pflip = pflip
        self.seed = seed
        self._rng = np.random.default_rng(seed)

    def load(
        self,
        name: str,
        source: Optional[str] = None,
        as_signed: bool = True,
        **kwargs,
    ) -> Union[SignedGraph, nx.Graph]:
        """
        Load a dataset by name.

        Parameters
        ----------
        name : str
            Dataset name.
        source : str, optional
            Data source. If None, tries builtin first, then searches.
        as_signed : bool, default True
            Whether to return a SignedGraph (True) or nx.Graph (False).
        **kwargs
            Additional arguments passed to specific loaders.

        Returns
        -------
        SignedGraph or nx.Graph
            Loaded network.
        """
        # Check built-in datasets first
        if name in BUILTIN_DATASETS or source == "builtin":
            G = self._load_builtin(name)
        elif source == "konect":
            G = self.load_konect(name, **kwargs)
        elif source == "snap":
            G = self.load_snap(name, **kwargs)
        elif source == "networkrepository":
            G = self.load_network_repository(name, **kwargs)
        elif source == "local" or source is None:
            # Check if it's a file path
            if Path(name).exists():
                G = self.load_from_file(name, **kwargs)
            elif name in BUILTIN_DATASETS:
                G = self._load_builtin(name)
            else:
                raise ValueError(
                    f"Unknown dataset '{name}'. Available built-in: "
                    f"{list(BUILTIN_DATASETS.keys())}"
                )
        else:
            raise ValueError(f"Unknown source '{source}'")

        if as_signed:
            return self._to_signed_graph(G)
        return G

    def _load_builtin(self, name: str) -> nx.Graph:
        """Load a built-in NetworkX dataset."""
        if name == "karate":
            return nx.karate_club_graph()
        elif name == "florentine_families":
            return nx.florentine_families_graph()
        elif name == "les_miserables":
            return nx.les_miserables_graph()
        elif name == "davis_southern_women":
            return nx.davis_southern_women_graph()
        elif name == "football":
            # Football is in GML format in NetworkX data
            try:
                url = (
                    "http://www-personal.umich.edu/~mejn/netdata/football.zip"
                )
                return self._load_gml_from_url(url, "football.gml")
            except Exception:
                # Fallback: create a similar structure
                warnings.warn(
                    "Could not download football network, using random substitute"
                )
                G = nx.random_partition_graph([12] * 10, 0.7, 0.05, seed=42)
                return G
        else:
            raise ValueError(f"Unknown built-in dataset: {name}")

    def _load_gml_from_url(self, url: str, filename: str) -> nx.Graph:
        """Load GML file from a URL (possibly zipped)."""
        import urllib.request
        import zipfile

        cache_path = self.cache_dir / filename
        if cache_path.exists():
            return nx.read_gml(cache_path)

        # Download
        with tempfile.TemporaryDirectory() as tmpdir:
            zip_path = Path(tmpdir) / "data.zip"
            urllib.request.urlretrieve(url, zip_path)

            with zipfile.ZipFile(zip_path, "r") as zf:
                # Extract GML file
                for name in zf.namelist():
                    if name.endswith(".gml"):
                        zf.extract(name, tmpdir)
                        gml_path = Path(tmpdir) / name
                        G = nx.read_gml(gml_path)
                        # Cache it
                        nx.write_gml(G, cache_path)
                        return G

        raise ValueError(f"No GML file found in {url}")

    def load_from_file(
        self,
        filepath: Union[str, Path],
        format: Optional[str] = None,
        delimiter: str = " ",
        comments: str = "#",
        nodetype: type = int,
        **kwargs,
    ) -> nx.Graph:
        """
        Load network from a local file.

        Parameters
        ----------
        filepath : str or Path
            Path to the file.
        format : str, optional
            File format. If None, inferred from extension.
            Supported: edgelist, adjlist, gml, graphml, pajek, gexf
        delimiter : str, default " "
            Delimiter for edge list files.
        comments : str, default "#"
            Comment character.
        nodetype : type, default int
            Type to convert node IDs to.
        **kwargs
            Additional arguments for the reader.

        Returns
        -------
        nx.Graph
            Loaded network.
        """
        filepath = Path(filepath)

        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")

        # Infer format from extension
        if format is None:
            suffix = filepath.suffix.lower()
            if suffix in (".gz", ".bz2"):
                # Handle compressed files
                suffix = Path(filepath.stem).suffix.lower()

            format_map = {
                ".edgelist": "edgelist",
                ".edges": "edgelist",
                ".txt": "edgelist",
                ".csv": "edgelist",
                ".adjlist": "adjlist",
                ".gml": "gml",
                ".graphml": "graphml",
                ".net": "pajek",
                ".gexf": "gexf",
                ".json": "json",
                ".mtx": "mtx",
            }
            format = format_map.get(suffix, "edgelist")

        # Load based on format
        if format == "edgelist":
            G = nx.read_edgelist(
                filepath,
                delimiter=delimiter,
                comments=comments,
                nodetype=nodetype,
                **kwargs,
            )
        elif format == "adjlist":
            G = nx.read_adjlist(filepath, comments=comments, nodetype=nodetype)
        elif format == "gml":
            G = nx.read_gml(filepath, **kwargs)
        elif format == "graphml":
            G = nx.read_graphml(filepath, **kwargs)
        elif format == "pajek":
            G = nx.read_pajek(filepath)
        elif format == "gexf":
            G = nx.read_gexf(filepath, **kwargs)
        elif format == "json":
            import json

            with open(filepath) as f:
                data = json.load(f)
            G = nx.node_link_graph(data)
        elif format == "mtx":
            G = self._read_mtx(filepath)
        else:
            raise ValueError(f"Unsupported format: {format}")

        return G

    def _read_mtx(self, filepath: Path) -> nx.Graph:
        """Read Matrix Market format."""
        try:
            from scipy.io import mmread

            A = mmread(filepath)
            G = nx.from_scipy_sparse_array(A)
            return G
        except ImportError:
            raise ImportError("scipy required for MTX format")

    def load_konect(
        self,
        name: str,
        **kwargs,
    ) -> nx.Graph:
        """
        Load dataset from KONECT (Koblenz Network Collection).

        Parameters
        ----------
        name : str
            Dataset name on KONECT (e.g., "dolphins", "karate").

        Returns
        -------
        nx.Graph
            Loaded network.

        Notes
        -----
        KONECT URL format: http://konect.cc/networks/{name}/
        Download URL: http://konect.cc/files/download.tsv.{name}.tar.bz2
        """
        import urllib.request
        import tarfile

        cache_path = self.cache_dir / f"konect_{name}"

        if cache_path.exists():
            # Load from cache
            edgelist = list(cache_path.glob("out.*"))
            if edgelist:
                return self._load_konect_format(edgelist[0])

        if not self.auto_download:
            raise ValueError(
                f"Dataset '{name}' not in cache. Set auto_download=True to download."
            )

        # Download from KONECT
        url = f"http://konect.cc/files/download.tsv.{name}.tar.bz2"

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                tar_path = Path(tmpdir) / f"{name}.tar.bz2"
                urllib.request.urlretrieve(url, tar_path)

                with tarfile.open(tar_path, "r:bz2") as tf:
                    tf.extractall(tmpdir)

                # Find the extracted directory
                extracted = [
                    p for p in Path(tmpdir).iterdir() if p.is_dir()
                ]
                if extracted:
                    # Move to cache
                    import shutil

                    shutil.move(str(extracted[0]), str(cache_path))

                    edgelist = list(cache_path.glob("out.*"))
                    if edgelist:
                        return self._load_konect_format(edgelist[0])

        except Exception as e:
            raise RuntimeError(f"Failed to download from KONECT: {e}")

        raise ValueError(f"Could not load KONECT dataset: {name}")

    def _load_konect_format(self, filepath: Path) -> nx.Graph:
        """Load KONECT edge list format."""
        G = nx.Graph()

        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("%"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    u, v = int(parts[0]), int(parts[1])
                    weight = float(parts[2]) if len(parts) > 2 else 1.0
                    G.add_edge(u, v, weight=weight)

        return G

    def load_snap(
        self,
        name: str,
        **kwargs,
    ) -> nx.Graph:
        """
        Load dataset from SNAP (Stanford Network Analysis Project).

        Parameters
        ----------
        name : str
            Dataset name on SNAP.

        Returns
        -------
        nx.Graph
            Loaded network.
        """
        import urllib.request

        cache_path = self.cache_dir / f"snap_{name}.txt.gz"

        # Common SNAP datasets
        snap_urls = {
            "facebook": "https://snap.stanford.edu/data/facebook_combined.txt.gz",
            "twitter": "https://snap.stanford.edu/data/twitter_combined.txt.gz",
            "epinions": "https://snap.stanford.edu/data/soc-Epinions1.txt.gz",
            "slashdot": "https://snap.stanford.edu/data/soc-Slashdot0811.txt.gz",
            "wiki-vote": "https://snap.stanford.edu/data/wiki-Vote.txt.gz",
            "email-enron": "https://snap.stanford.edu/data/email-Enron.txt.gz",
            "ca-grqc": "https://snap.stanford.edu/data/ca-GrQc.txt.gz",
            "ca-hepth": "https://snap.stanford.edu/data/ca-HepTh.txt.gz",
            "amazon": "https://snap.stanford.edu/data/amazon0302.txt.gz",
        }

        if name not in snap_urls:
            raise ValueError(
                f"Unknown SNAP dataset: {name}. "
                f"Available: {list(snap_urls.keys())}"
            )

        if not cache_path.exists():
            if not self.auto_download:
                raise ValueError(
                    f"Dataset '{name}' not in cache. "
                    "Set auto_download=True to download."
                )
            urllib.request.urlretrieve(snap_urls[name], cache_path)

        # Load gzipped edge list
        G = nx.Graph()
        with gzip.open(cache_path, "rt") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        u, v = int(parts[0]), int(parts[1])
                        G.add_edge(u, v)
                    except ValueError:
                        continue

        return G

    def load_network_repository(
        self,
        name: str,
        **kwargs,
    ) -> nx.Graph:
        """
        Load dataset from Network Repository.

        Parameters
        ----------
        name : str
            Dataset name.

        Returns
        -------
        nx.Graph
            Loaded network.
        """
        # Network Repository uses MTX format
        import urllib.request

        cache_path = self.cache_dir / f"nr_{name}.mtx"

        url = f"https://nrvis.com/download/data/soc/{name}.mtx"

        if not cache_path.exists():
            if not self.auto_download:
                raise ValueError(
                    f"Dataset '{name}' not in cache. "
                    "Set auto_download=True to download."
                )
            try:
                urllib.request.urlretrieve(url, cache_path)
            except Exception as e:
                raise RuntimeError(
                    f"Failed to download from Network Repository: {e}"
                )

        return self._read_mtx(cache_path)

    def _to_signed_graph(self, G: nx.Graph) -> SignedGraph:
        """
        Convert a NetworkX graph to SignedGraph.

        Assigns signs to edges and optionally flips some.
        """
        # Ensure integer node labels
        if not all(isinstance(n, int) for n in G.nodes()):
            G = nx.convert_node_labels_to_integers(G)

        # Add sign attributes
        for u, v in G.edges():
            if "sign" not in G[u][v]:
                G[u][v]["sign"] = self.default_sign
            if "weight" not in G[u][v]:
                G[u][v]["weight"] = G[u][v]["sign"]
            else:
                # Preserve weight sign if it exists
                if G[u][v]["weight"] < 0:
                    G[u][v]["sign"] = -1

        # Create SignedGraph
        sg = SignedGraph(G, pflip=self.pflip, seed=self.seed)

        # Apply random sign flipping if requested
        if self.pflip > 0:
            sg.flip_random_fract_edges()

        return sg

    def list_datasets(
        self,
        source: Optional[str] = None,
        tags: Optional[list[str]] = None,
    ) -> list[DatasetInfo]:
        """
        List available datasets.

        Parameters
        ----------
        source : str, optional
            Filter by source.
        tags : list of str, optional
            Filter by tags.

        Returns
        -------
        list of DatasetInfo
            Matching datasets.
        """
        datasets = list(BUILTIN_DATASETS.values())

        if source is not None:
            datasets = [d for d in datasets if d.source == source]

        if tags is not None:
            datasets = [
                d for d in datasets if any(t in d.tags for t in tags)
            ]

        return datasets

    def get_info(self, name: str) -> DatasetInfo:
        """
        Get metadata for a dataset.

        Parameters
        ----------
        name : str
            Dataset name.

        Returns
        -------
        DatasetInfo
            Dataset metadata.
        """
        if name in BUILTIN_DATASETS:
            return BUILTIN_DATASETS[name]
        raise ValueError(f"Unknown dataset: {name}")

    def clear_cache(self, name: Optional[str] = None) -> None:
        """
        Clear cached datasets.

        Parameters
        ----------
        name : str, optional
            Specific dataset to clear. If None, clears all.
        """
        if name is None:
            import shutil

            for path in self.cache_dir.iterdir():
                if path.is_dir():
                    shutil.rmtree(path)
                else:
                    path.unlink()
        else:
            # Find and remove specific dataset
            for path in self.cache_dir.iterdir():
                if name in path.name:
                    if path.is_dir():
                        import shutil

                        shutil.rmtree(path)
                    else:
                        path.unlink()
