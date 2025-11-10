import os
import sys
import tempfile
import types
import unittest
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
nx = pytest.importorskip("networkx")


class TestContactProcess(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._tmp_dir = tempfile.TemporaryDirectory()
        base_dir = Path(cls._tmp_dir.name)
        cls.data_dir = base_dir / "data"
        cls.log_dir = base_dir / "log"
        cls.data_dir.mkdir(parents=True, exist_ok=True)
        cls.log_dir.mkdir(parents=True, exist_ok=True)
        cls.ccore_dir = base_dir / "Ccore" / "bin"
        cls.ccore_dir.mkdir(parents=True, exist_ok=True)
        os.environ["LRGSG_CCORE_BIN"] = str(cls.ccore_dir)

        cls.src_path = Path(__file__).resolve().parents[1] / "src"
        sys.path.insert(0, str(cls.src_path))

        env_module = types.ModuleType("lrgsglib.config.lrgsg_env")
        env_module.LRGSG_LLIB = str(cls.src_path / "lrgsglib")
        env_module.LRGSG_DATA = str(cls.data_dir)
        env_module.LRGSG_LOG = str(cls.log_dir)
        env_module.LRGSG_CCORE_BIN = str(cls.ccore_dir)
        sys.modules["lrgsglib.config.lrgsg_env"] = env_module

        from lrgsglib.nx_patches import SignedGraph
        from lrgsglib.statsys.ContactProcess import ContactProcess

        cls.SignedGraph = SignedGraph
        cls.ContactProcess = ContactProcess

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp_dir.cleanup()
        if str(cls.src_path) in sys.path:
            sys.path.remove(str(cls.src_path))

    def _make_graph(self):
        G = nx.Graph()
        edges = [
            (0, 1, 1.0),
            (1, 2, -1.0),
            (2, 3, 1.0),
            (3, 0, 0.5),
        ]
        for u, v, w in edges:
            G.add_edge(u, v, weight=w)
        return self.SignedGraph(
            G,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )

    def test_init_contact_dynamics_uses_binary_state(self):
        sg = self._make_graph()
        model = self.ContactProcess(sg, runlang="py", ic="uniform", seed=123, simpref=1)
        model.init_contact_dynamics()
        try:
            unique = np.unique(model.s)
            self.assertTrue(np.all(np.isin(unique, (0, 1))))
        finally:
            model.reset_observables()

    def test_c1_builder_exports_and_arguments(self):
        sg = self._make_graph()
        model = self.ContactProcess(
            sg,
            runlang="C1",
            gamma=0.5,
            activation="relu",
            seed=0,
            simpref=1,
        )
        model.init_contact_dynamics()
        try:
            self.assertTrue(model.sfout.exists())
            self.assertTrue(hasattr(model.sg, "path_exp_edgl") and model.sg.path_exp_edgl.exists())
            self.assertTrue(hasattr(model.sg, "path_exp_adj") and model.sg.path_exp_adj.exists())
            expected_binary = self.ccore_dir / "ContactSimulator1"
            self.assertEqual(Path(model.cprogram[0]), expected_binary)
            self.assertEqual(model.cprogram[1], str(model.N))
            self.assertEqual(model.cprogram[2], f"{model.sg.pflip:.12g}")
            self.assertEqual(model.cprogram[3], f"{model.gamma:.12g}")
            self.assertEqual(model.cprogram[-1], model.activation)
        finally:
            if model.stderr_fopen and not model.stderr_fopen.closed:
                model.stderr_fopen.close()
            model.remove_run_c_files(remove_stderr=True)
            model.sg.remove_exported_files()


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
