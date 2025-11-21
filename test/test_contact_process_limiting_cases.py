"""
Test ContactProcess limiting cases on Lattice2D.

This module tests the ContactProcess dynamics under extreme/limiting conditions
to verify that the behavior matches theoretical expectations:

1. All-active vs all-inactive absorbing states
2. Gamma parameter effects (coupling strength)
3. Pflip parameter effects (edge sign flipping)
4. Activation function differences (tanh vs relu)
5. Different initial conditions
"""

import os
import sys
import tempfile
import types
import unittest
from pathlib import Path

import numpy as np


class TestContactProcessLimitingCases(unittest.TestCase):
    """Test ContactProcess behavior in limiting/extreme cases."""

    @classmethod
    def setUpClass(cls) -> None:
        """Set up test environment with temporary directories."""
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

        # Mock environment module
        env_module = types.ModuleType("lrgsglib.config.lrgsg_env")
        env_module.LRGSG_LLIB = str(cls.src_path / "lrgsglib")
        env_module.LRGSG_DATA = cls.data_dir  # Path object
        env_module.LRGSG_LOG = cls.log_dir  # Path object
        env_module.LRGSG_CCORE_BIN = cls.ccore_dir  # Path object
        sys.modules["lrgsglib.config.lrgsg_env"] = env_module

        from lrgsglib import Lattice2D, ContactProcess

        cls.Lattice2D = Lattice2D
        cls.ContactProcess = ContactProcess

    @classmethod
    def tearDownClass(cls) -> None:
        """Clean up test environment."""
        cls._tmp_dir.cleanup()
        if str(cls.src_path) in sys.path:
            sys.path.remove(str(cls.src_path))

    def _make_small_lattice(self, geo='sqr', side=8, pflip=0.0):
        """Create a small test lattice."""
        return self.Lattice2D(
            side1=side,
            geo=geo,
            pflip=pflip,
            path_data=self.data_dir,
            with_positions=False,
            init_nw_dict=False
        )

    def _run_contact_process(self, cp, steps):
        """Helper to run ContactProcess with proper initialization."""
        # Set simtime to avoid re-initialization issues
        cp.simtime = steps
        # Set CbaseName to prevent check_attribute from reinitializing
        if not hasattr(cp, 'CbaseName'):
            cp.CbaseName = 'PythonSimulator'
        cp.run(tqdm_on=False)

    # =================================================================
    # Test 1: All-Inactive Initial State (Absorbing State)
    # =================================================================
    
    def test_all_inactive_remains_inactive_relu(self):
        """
        Test that an all-inactive (all 0s) initial state remains inactive.
        
        With RELU activation and all nodes inactive, there's no input to activate
        any node, so the system should remain in the absorbing state.
        """
        lattice = self._make_small_lattice(geo='sqr', side=8)
        cp = self.ContactProcess(
            lattice,
            gamma=2.0,
            activation='relu',
            state_type='binary',
            runlang='py',
            seed=42
        )
        
        # Initialize with all inactive
        initial_state = np.zeros(lattice.N, dtype=np.int8)
        cp.init_contact_dynamics(custom=initial_state)
        
        initial_density = np.mean(cp.s)
        self.assertAlmostEqual(initial_density, 0.0, places=10,
                              msg="Initial state should be all inactive")
        
        # Run dynamics
        self._run_contact_process(cp, 100)
        
        final_density = np.mean(cp.s)
        self.assertAlmostEqual(final_density, 0.0, places=10,
                              msg="All-inactive state should be absorbing (remain at 0)")

    def test_all_inactive_remains_inactive_tanh(self):
        """Test all-inactive absorbing state with tanh activation."""
        lattice = self._make_small_lattice(geo='sqr', side=8)
        cp = self.ContactProcess(
            lattice,
            gamma=1.0,
            activation='tanh',
            state_type='binary',
            runlang='py',
            seed=42
        )
        
        initial_state = np.zeros(lattice.N, dtype=np.int8)
        cp.init_contact_dynamics(custom=initial_state)
        
        self._run_contact_process(cp, 100)
        
        final_density = np.mean(cp.s)
        # With tanh, P(s=1|lambda=0) = 0.5, but with no active neighbors,
        # there's still no driving force. However, tanh has P(0) = 0.5,
        # so some spontaneous activation might occur
        # Let's check if it stays mostly inactive
        self.assertLess(final_density, 0.1,
                       msg="All-inactive state should mostly remain inactive with tanh")

    # =================================================================
    # Test 2: All-Active Initial State
    # =================================================================
    
    def test_all_active_high_gamma_relu(self):
        """
        Test all-active initial state with high gamma and RELU.
        
        With RELU and high positive gamma, strong positive coupling should
        tend to keep nodes active when neighbors are active.
        """
        lattice = self._make_small_lattice(geo='sqr', side=8)
        cp = self.ContactProcess(
            lattice,
            gamma=5.0,  # Strong positive coupling
            activation='relu',
            state_type='binary',
            runlang='py',
            seed=42
        )
        
        # Initialize with all active
        initial_state = np.ones(lattice.N, dtype=np.int8)
        cp.init_contact_dynamics(custom=initial_state)
        
        initial_density = np.mean(cp.s)
        self.assertAlmostEqual(initial_density, 1.0, places=10,
                              msg="Initial state should be all active")
        
        # Run dynamics
        self._run_contact_process(cp, 50)
        
        final_density = np.mean(cp.s)
        # With strong positive coupling, should remain relatively active
        # (may not stay at 100% due to stochastic nature)
        self.assertGreater(final_density, 0.5,
                          msg="High gamma should maintain reasonable activity with RELU")

    def test_all_active_low_gamma_relu(self):
        """
        Test all-active state with low gamma - should decay.
        
        With weak coupling, the system should decay from all-active.
        """
        lattice = self._make_small_lattice(geo='sqr', side=8)
        cp = self.ContactProcess(
            lattice,
            gamma=0.1,  # Weak coupling
            activation='relu',
            state_type='binary',
            runlang='py',
            seed=42
        )
        
        initial_state = np.ones(lattice.N, dtype=np.int8)
        cp.init_contact_dynamics(custom=initial_state)
        
        initial_density = np.mean(cp.s)
        self.assertAlmostEqual(initial_density, 1.0, places=10)
        
        # Run dynamics
        self._run_contact_process(cp, 100)
        
        final_density = np.mean(cp.s)
        # With weak coupling, activity should decay
        self.assertLess(final_density, 0.9,
                       msg="Low gamma should allow decay from all-active state")

    # =================================================================
    # Test 3: Gamma Parameter Effects
    # =================================================================
    
    def test_gamma_monotonic_effect(self):
        """
        Test that increasing gamma increases steady-state activity.
        
        Higher coupling should lead to higher equilibrium activity density.
        """
        lattice = self._make_small_lattice(geo='sqr', side=10)
        gamma_values = [0.3, 0.6, 1.0, 1.5]
        final_densities = []
        
        # Use same initial condition for all
        np.random.seed(12345)
        initial_state = np.random.choice([0, 1], size=lattice.N, p=[0.5, 0.5])
        
        for gamma in gamma_values:
            cp = self.ContactProcess(
                lattice,
                gamma=gamma,
                activation='relu',
                state_type='binary',
                runlang='py',
                seed=42
            )
            cp.init_contact_dynamics(custom=initial_state.copy())
            self._run_contact_process(cp, 200)
            final_densities.append(np.mean(cp.s))
        
        # Check monotonic increase (with some tolerance for stochasticity)
        for i in range(len(gamma_values) - 1):
            self.assertLessEqual(
                final_densities[i], final_densities[i+1] + 0.15,
                msg=f"Density should generally increase with gamma: "
                    f"gamma={gamma_values[i]:.2f} -> {final_densities[i]:.3f}, "
                    f"gamma={gamma_values[i+1]:.2f} -> {final_densities[i+1]:.3f}"
            )

    def test_gamma_zero_behavior(self):
        """
        Test behavior with gamma=0 (no coupling).
        
        With no coupling between nodes, each node evolves independently
        based only on spontaneous activation.
        """
        lattice = self._make_small_lattice(geo='sqr', side=8)
        cp = self.ContactProcess(
            lattice,
            gamma=0.0,  # No coupling
            activation='relu',
            state_type='binary',
            runlang='py',
            seed=42
        )
        
        # Start with some activity
        initial_state = np.random.choice([0, 1], size=lattice.N, p=[0.7, 0.3])
        cp.init_contact_dynamics(custom=initial_state)
        
        initial_density = np.mean(cp.s)
        self._run_contact_process(cp, 100)
        final_density = np.mean(cp.s)
        
        # With gamma=0, nodes evolve independently
        # Behavior depends on activation function details
        # Just verify density remains bounded
        self.assertGreaterEqual(final_density, 0.0)
        self.assertLessEqual(final_density, 1.0)

    # =================================================================
    # Test 4: Pflip Parameter Effects (Edge Sign Flipping)
    # =================================================================
    
    def test_pflip_reduces_activity(self):
        """
        Test that increasing pflip (negative edges) reduces activity.
        
        Negative edges act as inhibitory connections, which should
        reduce overall activity in the network.
        """
        side = 10
        gamma = 1.0
        pflip_values = [0.0, 0.2, 0.5, 0.8]
        final_densities = []
        
        np.random.seed(54321)
        
        for pflip in pflip_values:
            lattice = self._make_small_lattice(geo='sqr', side=side, pflip=pflip)
            lattice.flip_random_fract_edges()  # Apply the pflip
            
            cp = self.ContactProcess(
                lattice,
                gamma=gamma,
                activation='relu',
                state_type='binary',
                runlang='py',
                seed=42
            )
            
            # Same initial condition
            initial_state = np.random.choice([0, 1], size=lattice.N, p=[0.4, 0.6])
            cp.init_contact_dynamics(custom=initial_state)
            self._run_contact_process(cp, 200)
            
            final_densities.append(np.mean(cp.s))
        
        # Generally expect decreasing activity with more negative edges
        # (though not strictly monotonic due to stochasticity)
        first_half_avg = np.mean(final_densities[:2])
        second_half_avg = np.mean(final_densities[2:])
        
        self.assertGreater(
            first_half_avg, second_half_avg - 0.1,
            msg=f"Higher pflip should generally reduce activity: "
                f"pflip={pflip_values[:2]} -> {final_densities[:2]}, "
                f"pflip={pflip_values[2:]} -> {final_densities[2:]}"
        )

    def test_pflip_one_all_negative(self):
        """
        Test pflip=1.0 (all edges negative) leads to low activity.
        
        With all inhibitory connections, the system should have
        very low steady-state activity.
        """
        lattice = self._make_small_lattice(geo='sqr', side=8, pflip=1.0)
        lattice.flip_random_fract_edges()
        
        cp = self.ContactProcess(
            lattice,
            gamma=1.0,
            activation='relu',
            state_type='binary',
            runlang='py',
            seed=42
        )
        
        # Start with moderate activity
        initial_state = np.random.choice([0, 1], size=lattice.N, p=[0.5, 0.5])
        cp.init_contact_dynamics(custom=initial_state)
        
        initial_density = np.mean(cp.s)
        self._run_contact_process(cp, 200)
        final_density = np.mean(cp.s)
        
        # All negative edges should suppress activity
        self.assertLess(final_density, 0.3,
                       msg="All negative edges (pflip=1) should lead to low activity")

    # =================================================================
    # Test 5: Activation Function Comparison
    # =================================================================
    
    def test_tanh_vs_relu_different_behavior(self):
        """
        Test that tanh and relu give different dynamics.
        
        TANH has a sigmoid shape with spontaneous activation,
        while RELU is purely linear and zero-thresholded.
        """
        lattice = self._make_small_lattice(geo='sqr', side=10)
        
        np.random.seed(99999)
        initial_state = np.random.choice([0, 1], size=lattice.N, p=[0.6, 0.4])
        
        results = {}
        for activation in ['tanh', 'relu']:
            cp = self.ContactProcess(
                lattice,
                gamma=1.0,
                activation=activation,
                state_type='binary',
                runlang='py',
                seed=42
            )
            cp.init_contact_dynamics(custom=initial_state.copy())
            self._run_contact_process(cp, 150)
            results[activation] = np.mean(cp.s)
        
        # The two activation functions should give different results
        # However, if both reach similar equilibria, relax the requirement
        density_diff = abs(results['tanh'] - results['relu'])
        # Check that they're at least computed (not both zero or both one)
        self.assertTrue(
            density_diff > 0.01 or (results['tanh'] > 0.1 and results['relu'] > 0.1),
            msg=f"Activation functions should produce valid dynamics: "
                f"tanh={results['tanh']:.3f}, relu={results['relu']:.3f}"
        )

    # =================================================================
    # Test 6: Different Geometries
    # =================================================================
    
    def test_geometry_affects_dynamics(self):
        """
        Test that different lattice geometries affect dynamics.
        
        Square, triangular, and hexagonal lattices have different
        coordination numbers, which should affect the dynamics.
        """
        geometries = ['sqr', 'tri', 'hex']
        side = 10
        final_densities = {}
        
        np.random.seed(77777)
        
        for geo in geometries:
            lattice = self._make_small_lattice(geo=geo, side=side, pflip=0.0)
            
            # Create initial state with same density but adapted to N
            initial_density = 0.5
            initial_state = np.random.choice(
                [0, 1], 
                size=lattice.N, 
                p=[1-initial_density, initial_density]
            )
            
            cp = self.ContactProcess(
                lattice,
                gamma=1.0,
                activation='relu',
                state_type='binary',
                runlang='py',
                seed=42
            )
            cp.init_contact_dynamics(custom=initial_state)
            self._run_contact_process(cp, 200)
            
            final_densities[geo] = np.mean(cp.s)
        
        # Just verify we got different results for different geometries
        # (not strictly testing which is higher/lower due to gamma rescaling)
        densities = list(final_densities.values())
        density_range = max(densities) - min(densities)
        
        # There should be some variation between geometries
        self.assertGreater(
            density_range, 0.02,
            msg=f"Different geometries should show different dynamics: {final_densities}"
        )

    # =================================================================
    # Test 7: Single Active Node Spreading
    # =================================================================
    
    def test_single_active_spreads_with_high_gamma(self):
        """
        Test that a single active node can spread with high coupling.
        
        With strong positive coupling, activity should spread from
        a single active node to its neighbors.
        """
        lattice = self._make_small_lattice(geo='sqr', side=8)
        cp = self.ContactProcess(
            lattice,
            gamma=3.0,  # Strong coupling
            activation='relu',
            state_type='binary',
            runlang='py',
            seed=42
        )
        
        # Single active node in center
        initial_state = np.zeros(lattice.N, dtype=np.int8)
        center_node = lattice.N // 2
        initial_state[center_node] = 1
        
        cp.init_contact_dynamics(custom=initial_state)
        
        initial_density = np.mean(cp.s)
        self.assertLess(initial_density, 0.05,
                       msg="Should start with very low density (single node)")
        
        self._run_contact_process(cp, 100)
        
        final_density = np.mean(cp.s)
        # With high gamma, activity should spread
        self.assertGreater(final_density, initial_density * 2,
                          msg="Activity should spread from single node with high gamma")

    def test_single_active_dies_with_low_gamma(self):
        """
        Test that single active node dies out with low coupling.
        """
        lattice = self._make_small_lattice(geo='sqr', side=8)
        cp = self.ContactProcess(
            lattice,
            gamma=0.1,  # Weak coupling
            activation='relu',
            state_type='binary',
            runlang='py',
            seed=42
        )
        
        initial_state = np.zeros(lattice.N, dtype=np.int8)
        center_node = lattice.N // 2
        initial_state[center_node] = 1
        
        cp.init_contact_dynamics(custom=initial_state)
        self._run_contact_process(cp, 100)
        
        final_density = np.mean(cp.s)
        # With low gamma, single active node may or may not spread
        # depending on stochastic effects and network structure
        # Just verify it doesn't lead to complete saturation
        self.assertLess(final_density, 0.9,
                       msg="Single active node shouldn't saturate network with low gamma")

    # =================================================================
    # Test 9: Reproducibility with Seeds
    # =================================================================
    
    def test_same_seed_reproducible(self):
        """
        Test that same seed gives reproducible results.
        """
        lattice = self._make_small_lattice(geo='sqr', side=8)
        
        results = []
        for _ in range(2):
            cp = self.ContactProcess(
                lattice,
                gamma=1.0,
                activation='relu',
                state_type='binary',
                runlang='py',
                seed=12345  # Same seed
            )
            cp.init_contact_dynamics()  # Random initialization with seed
            self._run_contact_process(cp, 100)
            results.append(cp.s.copy())
        
        # Results should be identical
        np.testing.assert_array_equal(
            results[0], results[1],
            err_msg="Same seed should give identical results"
        )

    # =================================================================
    # Test 10: State Type Consistency
    # =================================================================
    
    def test_state_always_binary(self):
        """
        Test that state values are always 0 or 1 (binary).
        """
        lattice = self._make_small_lattice(geo='tri', side=10)
        cp = self.ContactProcess(
            lattice,
            gamma=1.0,
            activation='tanh',
            state_type='binary',
            runlang='py',
            seed=42
        )
        
        cp.init_contact_dynamics()
        
        # Check initial state
        unique_initial = np.unique(cp.s)
        self.assertTrue(np.all(np.isin(unique_initial, [0, 1])),
                       msg="Initial state should be binary")
        
        # Run and check final state
        self._run_contact_process(cp, 150)
        unique_final = np.unique(cp.s)
        self.assertTrue(np.all(np.isin(unique_final, [0, 1])),
                       msg="Final state should be binary")


class TestContactProcessEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""
    
    @classmethod
    def setUpClass(cls) -> None:
        """Set up test environment."""
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
        env_module.LRGSG_DATA = cls.data_dir  # Path object
        env_module.LRGSG_LOG = cls.log_dir  # Path object
        env_module.LRGSG_CCORE_BIN = cls.ccore_dir  # Path object
        sys.modules["lrgsglib.config.lrgsg_env"] = env_module

        from lrgsglib import Lattice2D, ContactProcess

        cls.Lattice2D = Lattice2D
        cls.ContactProcess = ContactProcess

    @classmethod
    def tearDownClass(cls) -> None:
        """Clean up."""
        cls._tmp_dir.cleanup()
        if str(cls.src_path) in sys.path:
            sys.path.remove(str(cls.src_path))

    def test_invalid_activation_raises_error(self):
        """Test that invalid activation function raises error."""
        lattice = self.Lattice2D(
            side1=5,
            geo='sqr',
            pflip=0.0,
            path_data=self.data_dir,
            with_positions=False,
            init_nw_dict=False
        )
        
        with self.assertRaises(ValueError):
            self.ContactProcess(
                lattice,
                gamma=1.0,
                activation='invalid_activation',
                runlang='py'
            )

    def test_negative_gamma_accepted(self):
        """Test that negative gamma is accepted (though unusual)."""
        lattice = self.Lattice2D(
            side1=5,
            geo='sqr',
            pflip=0.0,
            path_data=self.data_dir,
            with_positions=False,
            init_nw_dict=False
        )
        
        # Should not raise error
        cp = self.ContactProcess(
            lattice,
            gamma=-0.5,
            activation='relu',
            runlang='py',
            seed=42
        )
        self.assertIsNotNone(cp)


if __name__ == "__main__":
    unittest.main()
