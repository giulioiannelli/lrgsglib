import os
import sys
import types
from pathlib import Path

import pytest


@pytest.mark.filterwarnings("ignore:Animation was deleted without rendering")
def test_collect_frames_and_make_animation(tmp_path):
    mpl = pytest.importorskip("matplotlib")
    mpl.use("Agg")
    pytest.importorskip("numpy")
    pytest.importorskip("networkx")

    data_dir = tmp_path / "data"
    log_dir = tmp_path / "log"
    ccore_dir = tmp_path / "Ccore" / "bin"
    data_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    ccore_dir.mkdir(parents=True, exist_ok=True)
    os.environ["LRGSG_CCORE_BIN"] = str(ccore_dir)

    env_module = types.ModuleType("lrgsglib.config.lrgsg_env")
    env_module.LRGSG_LLIB = ""
    env_module.LRGSG_DATA = str(data_dir)
    env_module.LRGSG_LOG = log_dir
    env_module.LRGSG_CCORE_BIN = str(ccore_dir)
    sys.modules["lrgsglib.config.lrgsg_env"] = env_module

    from lrgsglib.animations import collect_contact_process_frames
    from lrgsglib.graphs.nx import Lattice2DNX as Lattice2D
    from lrgsglib.statsys.ContactProcess import ContactProcessSIR

    lattice = Lattice2D(
        side1=4,
        side2=4,
        geo="squared",
        pbc=False,
        init_nw_dict=False,
        path_data=data_dir,
        path_plot=log_dir,
    )

    cp = ContactProcessSIR(
        lattice,
        mu=0.2,
        runlang="py",
        ic="delta",
        state_type="binary",
        savedyn=False,
        seed=0,
        steps=1,
    )
    cp.init_contact_dynamics()

    out = collect_contact_process_frames(cp, steps=3, backend="py")
    assert len(out.frames) == 4  # t0 + 3 sweeps

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    result = lattice.make_animation(fig, ax, out.frames, add_colorbar=False)
    assert hasattr(result.animation, "save")
    plt.close(fig)
