import os
import sys
import subprocess
import tempfile
import types
import unittest
from pathlib import Path

import pytest

nx = pytest.importorskip("networkx")
np = pytest.importorskip("numpy")
from unittest import mock

class TestContactProcess(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._tmp_dir = tempfile.TemporaryDirectory()
        base_dir = Path(cls._tmp_dir.name)
        cls.data_dir = base_dir / "data"
        cls.log_dir = base_dir / "log"
        cls.data_dir.mkdir(parents=True, exist_ok=True)
        cls.log_dir.mkdir(parents=True, exist_ok=True)
        # sys.path needed for kernels/parsers imports (not installed packages)
        cls.src_path = Path(__file__).resolve().parents[2] / "src"
        sys.path.insert(0, str(cls.src_path))

        env_module = types.ModuleType("lrgsglib.config.lrgsg_env")
        env_module.LRGSG_LLIB = ""
        env_module.LRGSG_DATA = str(cls.data_dir)
        env_module.LRGSG_LOG = str(cls.log_dir)
        sys.modules["lrgsglib.config.lrgsg_env"] = env_module

        from lrgsglib.graphs.nx import SignedGraphNX as SignedGraph
        from lrgsglib.statsys.ContactProcess import ContactProcessEI, ContactProcessSIR
        from kernels.ContactProcessDynamics import initialize_contact_process_dict_args
        from kernels import L2D_ContactProcess as cp_kernel
        from parsers import L2D_ContactProcess as cp_parser

        cls.SignedGraph = SignedGraph
        cls.ContactProcessSIR = ContactProcessSIR
        cls.ContactProcessEI = ContactProcessEI
        cls.initialize_contact_process_dict_args = staticmethod(
            initialize_contact_process_dict_args
        )
        cls.cp_kernel = cp_kernel
        cls.cp_parser = cp_parser

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
        model = self.ContactProcessSIR(
            sg,
            runlang="py",
            ic="uniform",
            seed=123,
            steps=1,
        )
        model.init_contact_dynamics()
        try:
            unique = np.unique(model.s)
            self.assertTrue(np.all(np.isin(unique, (0, 1))))
        finally:
            model.reset_observables()

    def test_c1_builder_exports_and_arguments(self):
        path_graph = nx.path_graph(3)
        for u, v in path_graph.edges:
            path_graph[u][v]["weight"] = 1.0
        sg = self.SignedGraph(
            path_graph,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        model = self.ContactProcessEI(
            sg,
            runlang="C1c",
            gamma=0.6,
            activation="relu",
            num_log_samples=25,
            seed=0,
            steps=1,
        )
        model._set_time_controls(steps=42)
        model.init_contact_dynamics()
        try:
            self.assertTrue(model.sfout.exists())
            self.assertTrue(
                hasattr(model.sg, "path_exp_edgl")
                and model.sg.path_exp_edgl.exists()
            )
            self.assertEqual(Path(model.cprogram[0]).name, "ContactSimulator1c")
            self.assertEqual(model.cprogram[1], str(model.N))
            self.assertEqual(model.cprogram[2], f"{model.sg.pflip:.12g}")
            self.assertNotEqual(model.gamma_eff, model.gamma)
            self.assertEqual(model.cprogram[3], f"{model.gamma_eff:.12g}")
            self.assertEqual(model.cprogram[4], f"{model.steps}")
            self.assertEqual(model.cprogram[8], model.out_id)
            self.assertEqual(model.cprogram[-2], model.activation)
            self.assertEqual(model.cprogram[-1], f"{model.num_log_samples}")
        finally:
            if model.stderr_fopen and not model.stderr_fopen.closed:
                model.stderr_fopen.close()
            model.remove_run_c_files(remove_stderr=True)
            model.sg.remove_exported_files()

    def test_c1e_builder_arguments_include_cached_lambda_variant(self):
        path_graph = nx.path_graph(4)
        for u, v in path_graph.edges:
            path_graph[u][v]["weight"] = 1.0
        sg = self.SignedGraph(
            path_graph,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        model = self.ContactProcessEI(
            sg,
            runlang="C1e",
            gamma=0.5,
            activation="tanh",
            num_log_samples=15,
            seed=0,
            steps=1,
        )
        model._set_time_controls(steps=30)
        model.init_contact_dynamics()
        try:
            self.assertEqual(Path(model.cprogram[0]).name, "ContactSimulator1e")
            self.assertEqual(model.cprogram[1], str(model.N))
            self.assertEqual(model.cprogram[2], f"{model.sg.pflip:.12g}")
            self.assertEqual(model.cprogram[3], f"{model.gamma_eff:.12g}")
            self.assertEqual(model.cprogram[4], f"{model.steps}")
            self.assertEqual(model.cprogram[-2], model.activation)
            self.assertEqual(model.cprogram[-1], f"{model.num_log_samples}")
        finally:
            if model.stderr_fopen and not model.stderr_fopen.closed:
                model.stderr_fopen.close()
            model.remove_run_c_files(remove_stderr=True)
            model.sg.remove_exported_files()

    def test_c1f_builder_arguments_include_gillespie_variant(self):
        path_graph = nx.path_graph(3)
        for u, v in path_graph.edges:
            path_graph[u][v]["weight"] = 1.0
        sg = self.SignedGraph(
            path_graph,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        model = self.ContactProcessEI(
            sg,
            runlang="C1f",
            gamma=0.75,
            activation="relu",
            num_log_samples=12,
            seed=1,
            steps=1,
        )
        model._set_time_controls(steps=25)
        model.init_contact_dynamics()
        try:
            self.assertEqual(Path(model.cprogram[0]).name, "ContactSimulator1f")
            self.assertEqual(model.cprogram[1], str(model.N))
            self.assertEqual(model.cprogram[2], f"{model.sg.pflip:.12g}")
            self.assertEqual(model.cprogram[3], f"{model.gamma_eff:.12g}")
            self.assertEqual(model.cprogram[4], f"{model.steps}")
            self.assertEqual(model.cprogram[-2], model.activation)
            self.assertEqual(model.cprogram[-1], f"{model.num_log_samples}")
        finally:
            if model.stderr_fopen and not model.stderr_fopen.closed:
                model.stderr_fopen.close()
            model.remove_run_c_files(remove_stderr=True)
            model.sg.remove_exported_files()

    def test_sir_step_recovery_and_infection(self):
        sg = self._make_graph()
        model = self.ContactProcessSIR(
            sg,
            runlang="py",
            mu=5.0,
            steps=1,
            seed=0,
        )

        # Force deterministic random draws
        random_calls = mock.Mock(side_effect=[0.0, 0.0])
        with mock.patch("numpy.random.random", random_calls):
            model.s = np.array([1, 0, 0, 0], dtype=np.int8)
            model.ds1step(0)  # high mu => deactivate
            self.assertEqual(model.s[0], 0)

            model.s = np.array([1, 0, 0, 0], dtype=np.int8)
            # Activate node 1 due to strong positive neighbour
            sg.gr[sg.on_g][0][1]["weight"] = 10.0
            model.ds1step(1)
            self.assertEqual(model.s[1], 1)

    def test_ei_run_allows_python_backend(self):
        sg = self._make_graph()
        model = self.ContactProcessEI(sg, gamma=0.4, runlang="py", steps=1, seed=0)
        model.init_contact_dynamics()
        model.run(tqdm_on=False, steps=2)
        unique = np.unique(model.s)
        self.assertTrue(np.all(np.isin(unique, (0, 1))))

    def test_c_backend_run_reads_stdout_state(self):
        sg = self._make_graph()
        model = self.ContactProcessSIR(sg, mu=0.2, runlang="C0", steps=1, seed=0)
        model._set_time_controls(steps=5)

        fake_binary = self.ccore_dir / "ContactSimulator0"
        fake_binary.write_text("#!/bin/sh\n")

        stdout_state = np.ones(model.N, dtype=np.int8)

        with mock.patch("subprocess.run") as run_mock:
            run_mock.return_value = subprocess.CompletedProcess(
                args=[str(fake_binary)], returncode=0, stdout=stdout_state.tobytes()
            )
            model.init_contact_dynamics()
            model.run(verbose=False, clean_export=False)

        np.testing.assert_array_equal(model.s, stdout_state)
        model.remove_run_c_files(remove_stderr=True)
        model.sg.remove_exported_files()

    def test_cli_validation_blocks_missing_gamma(self):
        with self.assertRaises(SystemExit):
            self.cp_parser.parse_arguments(
                self.cp_parser.parser,
                ["4", "0.1", "--dynamics", "EI", "--runlang", "C1c"],
            )

    def test_cli_validation_blocks_sir_c1(self):
        with self.assertRaises(SystemExit):
            self.cp_parser.parse_arguments(
                self.cp_parser.parser,
                ["4", "0.1", "--dynamics", "SIR", "--runlang", "C1c", "--mu", "0.5"],
            )

    def test_kernel_validates_backend_choice(self):
        args = types.SimpleNamespace(
            dynamics="EI",
            gamma=None,
            mu=None,
            activation="tanh",
            num_log_samples=10,
            state_type="binary",
            init_cond="uniform",
            runlang="py",
            steps=1,
            simref=None,
            randstr=False,
            out_suffix="",
        )
        with self.assertRaises(ValueError):
            self.initialize_contact_process_dict_args(args, "suffix")

    def test_l2d_contact_process_rejects_unwired_backends(self):
        args = self.cp_parser.parse_arguments(
            self.cp_parser.parser,
            ["4", "0.1", "--dynamics", "EI", "--runlang", "C1a", "--gamma", "0.3"],
        )
        with self.assertRaises(NotImplementedError):
            self.cp_kernel.run_simulation(args)

    def test_simref_scales_steps(self):
        sg = self._make_graph()
        model = self.ContactProcessEI(sg, gamma=0.5, simref=5)
        self.assertEqual(model.simref, 5)
        self.assertEqual(model.steps, 5 * model.N)

    def test_out_id_includes_gamma_and_random_suffix(self):
        sg = self._make_graph()
        model = self.ContactProcessEI(
            sg,
            gamma=0.499,
            runlang="C1c",
            simref=2,
            rndStr=True,
            out_suffix="uniform_rand",
        )
        model.init_contact_dynamics()
        try:
            self.assertTrue(model.run_id)
            self.assertIn("gamma=0.499", model.out_id)
            self.assertTrue(model.out_id.endswith(model.run_id))
        finally:
            if model.stderr_fopen and not model.stderr_fopen.closed:
                model.stderr_fopen.close()
            model.remove_run_c_files(remove_stderr=True)
            model.sg.remove_exported_files()

    def test_cu_unified_builder_arguments(self):
        """Test that CU runlang builds correct unified ContactSimulator args."""
        path_graph = nx.path_graph(4)
        for u, v in path_graph.edges:
            path_graph[u][v]["weight"] = 1.0
        sg = self.SignedGraph(
            path_graph,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        model = self.ContactProcessEI(
            sg,
            runlang="CU",
            gamma=0.5,
            activation="tanh",
            update_mode="cached",
            output_mode="density",
            num_log_samples=50,
            seed=0,
            steps=1,
        )
        model._set_time_controls(steps=100)
        model.init_contact_dynamics()
        try:
            # Unified binary: ContactSimulator (no suffix)
            self.assertEqual(Path(model.cprogram[0]).name, "ContactSimulator")
            self.assertEqual(model.cprogram[1], str(model.N))
            self.assertEqual(model.cprogram[2], f"{model.sg.pflip:.12g}")
            self.assertEqual(model.cprogram[3], f"{model.gamma_eff:.12g}")
            self.assertEqual(model.cprogram[4], f"{model.steps}")
            # datdir, syshape, run_id, out_id
            self.assertEqual(model.cprogram[8], model.out_id)
            # activation, update_mode, output_mode, num_samples
            self.assertEqual(model.cprogram[9], "tanh")
            self.assertEqual(model.cprogram[10], "cached")
            self.assertEqual(model.cprogram[11], "density")
            self.assertEqual(model.cprogram[12], "50")
        finally:
            if model.stderr_fopen and not model.stderr_fopen.closed:
                model.stderr_fopen.close()
            model.remove_run_c_files(remove_stderr=True)
            model.sg.remove_exported_files()

    def test_cu_unified_builder_final_mode_no_samples(self):
        """Test that CU with final output mode omits num_samples."""
        path_graph = nx.path_graph(3)
        for u, v in path_graph.edges:
            path_graph[u][v]["weight"] = 1.0
        sg = self.SignedGraph(
            path_graph,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        model = self.ContactProcessEI(
            sg,
            runlang="CU",
            gamma=0.3,
            activation="relu",
            update_mode="naive",
            output_mode="final",
            seed=1,
            steps=1,
        )
        model._set_time_controls(steps=50)
        model.init_contact_dynamics()
        try:
            self.assertEqual(Path(model.cprogram[0]).name, "ContactSimulator")
            self.assertEqual(model.cprogram[9], "relu")
            self.assertEqual(model.cprogram[10], "naive")
            self.assertEqual(model.cprogram[11], "final")
            # No num_samples for final mode
            self.assertEqual(len(model.cprogram), 12)
        finally:
            if model.stderr_fopen and not model.stderr_fopen.closed:
                model.stderr_fopen.close()
            model.remove_run_c_files(remove_stderr=True)
            model.sg.remove_exported_files()

    def test_update_mode_validation(self):
        """Test that invalid update_mode raises ValueError."""
        sg = self._make_graph()
        with self.assertRaises(ValueError):
            self.ContactProcessEI(
                sg,
                gamma=0.5,
                update_mode="invalid",
            )

    def test_output_mode_validation(self):
        """Test that invalid output_mode raises ValueError."""
        sg = self._make_graph()
        with self.assertRaises(ValueError):
            self.ContactProcessEI(
                sg,
                gamma=0.5,
                output_mode="invalid",
            )


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
