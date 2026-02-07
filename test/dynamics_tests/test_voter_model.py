import os
import sys
import tempfile
import types
import unittest
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
nx = pytest.importorskip("networkx")


class TestVoterModel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp_dir = tempfile.TemporaryDirectory()
        base_dir = Path(cls._tmp_dir.name)
        cls.data_dir = base_dir / "data"
        cls.log_dir = base_dir / "log"
        cls.data_dir.mkdir(parents=True, exist_ok=True)
        cls.log_dir.mkdir(parents=True, exist_ok=True)
        cls.ccore_dir = base_dir / "Ccore" / "bin"
        cls.ccore_dir.mkdir(parents=True, exist_ok=True)
        os.environ["LRGSG_CCORE_BIN"] = str(cls.ccore_dir)

        env_module = types.ModuleType("lrgsglib.config.lrgsg_env")
        env_module.LRGSG_LLIB = ""
        env_module.LRGSG_DATA = str(cls.data_dir)
        env_module.LRGSG_LOG = str(cls.log_dir)
        env_module.LRGSG_CCORE_BIN = str(cls.ccore_dir)
        sys.modules["lrgsglib.config.lrgsg_env"] = env_module

        from lrgsglib.graphs.nx import SignedGraphNX as SignedGraph
        from lrgsglib.statsys.VoterModel import VoterModel

        cls.SignedGraph = SignedGraph
        cls.VoterModel = VoterModel

    @classmethod
    def tearDownClass(cls):
        cls._tmp_dir.cleanup()

    def _make_graph(self):
        G = nx.Graph()
        edges = [
            (0, 1, 1.0),
            (1, 2, -1.0),
            (2, 3, 1.0),
            (3, 0, 1.0),
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

    def test_python_dynamics_collects_history(self):
        sg = self._make_graph()
        model = self.VoterModel(
            sg,
            runlang="py",
            savedyn=True,
            save_magnetization=True,
            eqSTEP=5,
            seed=123,
        )
        model.run(tqdm_on=False)
        self.assertEqual(len(model.magn), 5)
        self.assertEqual(len(model.s_t), 5)
        self.assertTrue(np.all(np.isin(model.s, (-1, 1))))

    def test_init_voter_dynamics_exports_required_files(self):
        sg = self._make_graph()
        model = self.VoterModel(sg, runlang="C", eqSTEP=2, seed=0)
        model.init_voter_dynamics()
        try:
            self.assertTrue(model.sfout.exists())
            self.assertTrue(hasattr(model.sg, "path_exp_edgl") and model.sg.path_exp_edgl.exists())
            self.assertTrue(hasattr(model.sg, "path_exp_adj") and model.sg.path_exp_adj.exists())
        finally:
            if model.stderr_fopen and not model.stderr_fopen.closed:
                model.stderr_fopen.close()
            model.remove_run_c_files(remove_stderr=True)
            model.sg.remove_exported_files()

    def test_run_c_raises_when_binary_missing(self):
        sg = self._make_graph()
        model = self.VoterModel(sg, runlang="C", eqSTEP=1, seed=0)
        model.init_voter_dynamics()
        with self.assertRaises(FileNotFoundError):
            model.run(tqdm_on=False, clean_export=False)
        if model.stderr_fopen and not model.stderr_fopen.closed:
            model.stderr_fopen.close()
        model.remove_run_c_files(remove_stderr=True)
        model.sg.remove_exported_files()


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
