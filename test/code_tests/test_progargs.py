"""Tests for config/progargs argument definitions.

Validates that argument dictionaries are well-formed and consistent.
"""
import pytest


class TestProgargsDictStructure:
    """Verify progargs dictionaries have valid argparse-compatible structure."""

    def _check_args_dict(self, args_dict):
        """Verify a positional args dict has valid structure."""
        for name, config in args_dict.items():
            assert isinstance(name, str), f"Arg name must be str, got {type(name)}"
            assert isinstance(config, dict), f"Arg config for {name} must be dict"
            assert "type" in config or "action" in config, (
                f"Arg {name} must have 'type' or 'action'"
            )

    def _check_opt_args_dict(self, opt_args_dict):
        """Verify an optional args dict has valid structure."""
        for flags, config in opt_args_dict.items():
            assert isinstance(flags, tuple), (
                f"Optional arg flags must be tuple, got {type(flags)}"
            )
            assert len(flags) >= 1, "Must have at least one flag"
            for flag in flags:
                assert isinstance(flag, str), f"Flag must be str, got {type(flag)}"
                assert flag.startswith("-"), f"Flag '{flag}' must start with '-'"
            assert isinstance(config, dict), f"Config for {flags} must be dict"

    def test_lattice2d_args(self):
        from lrgsglib.config.progargs.Lattice2D import L2D_args, L2D_opt_args

        self._check_args_dict(L2D_args)
        self._check_opt_args_dict(L2D_opt_args)

    def test_lattice3d_args(self):
        from lrgsglib.config.progargs.Lattice3D import L3D_args, L3D_opt_args

        self._check_args_dict(L3D_args)
        self._check_opt_args_dict(L3D_opt_args)

    def test_erdosrenyi_args(self):
        from lrgsglib.config.progargs.ErdosRenyi import ER_args, ER_opt_args

        self._check_args_dict(ER_args)
        self._check_opt_args_dict(ER_opt_args)

    def test_ising_args(self):
        from lrgsglib.config.progargs.IsingDynamics import IsDyn_args, IsDyn_opt_args

        self._check_args_dict(IsDyn_args)
        self._check_opt_args_dict(IsDyn_opt_args)

    def test_contact_process_args(self):
        from lrgsglib.config.progargs.ContactProcess import (
            Contact_args,
            Contact_opt_args,
        )

        # Contact_args uses tuple keys (optional-style flags)
        self._check_opt_args_dict(Contact_args)
        self._check_opt_args_dict(Contact_opt_args)

    def test_slaplspect_args(self):
        from lrgsglib.config.progargs.SlaplSpect import (
            SlaplSpect_common_optional_args_dict,
        )

        self._check_opt_args_dict(SlaplSpect_common_optional_args_dict)


class TestProgargDefaults:
    """Verify default values are sensible types."""

    def test_lattice2d_defaults(self):
        from lrgsglib.config.progargs.Lattice2D import L2D_opt_args

        for flags, config in L2D_opt_args.items():
            if "default" in config and "type" in config:
                default = config["default"]
                expected_type = config["type"]
                if default is not None:
                    assert isinstance(default, expected_type), (
                        f"Default for {flags} is {type(default)}, "
                        f"expected {expected_type}"
                    )

    def test_ising_defaults(self):
        from lrgsglib.config.progargs.IsingDynamics import IsDyn_opt_args as Ising_opt_args

        for flags, config in Ising_opt_args.items():
            if "default" in config and "type" in config:
                default = config["default"]
                expected_type = config["type"]
                if default is not None:
                    assert isinstance(default, expected_type), (
                        f"Default for {flags} is {type(default)}, "
                        f"expected {expected_type}"
                    )
