from .common import *
from .BinDynSys import BinDynSys


class ContactProcess(BinDynSys):
    """Simple contact process dynamics on a weighted graph."""

    def __init__(self, sg: SignedGraph, mu: float = 1.0, **kwargs) -> None:
        super().__init__(sg, **kwargs)
        self.mu = mu
        self.dynpath = getattr(self.sg, "path_lrgsg", self.sg.path_data)

    def init_state(self, custom=None):
        match self.ic:
            case "uniform" | "random" | "rand":
                self.s = np.random.choice([0, 1], size=self.sg.N).astype(
                    np.int8
                )
            case "homogeneous" | "homo":
                self.s = np.ones(self.sg.N, dtype=np.int8)
            case "delta":
                self.s = np.zeros(self.sg.N, dtype=np.int8)
                self.s[np.random.randint(self.sg.N)] = 1
            case "custom":
                self.s = custom.astype(np.int8)
            case _:
                raise ValueError("Invalid initial condition.")

    def infection_rate(self, j: int, on_g: str = SG_GRAPH_REPR) -> float:
        neigh = (
            self.sg.gr[on_g].predecessors(j)
            if hasattr(self.sg.gr[on_g], "predecessors")
            else self.sg.gr[on_g].neighbors(j)
        )
        rate = 0.0
        for i in neigh:
            w = self.sg.gr[on_g][i][j].get("weight", 1.0)
            rate += self.s[i] * max(w, 0.0)
        return rate

    def recovery_rate(self, j: int, on_g: str = SG_GRAPH_REPR) -> float:
        neigh = (
            self.sg.gr[on_g].predecessors(j)
            if hasattr(self.sg.gr[on_g], "predecessors")
            else self.sg.gr[on_g].neighbors(j)
        )
        rate = self.mu
        for i in neigh:
            w = self.sg.gr[on_g][i][j].get("weight", 1.0)
            rate += self.s[i] * max(-w, 0.0)
        return rate

    def ds1step(self, j: int):
        if self.s[j]:
            r = self.recovery_rate(j)
            if np.random.rand() < 1 - np.exp(-r):
                self.s[j] = 0
        else:
            r = self.infection_rate(j)
            if np.random.rand() < 1 - np.exp(-r):
                self.s[j] = 1

    def run_py(self):
        dsNstep = self.dsNstep()
        nodes = list(self.sg.gr[self.sg.on_g].nodes())
        for _ in range(self.simtime):
            smp = random.sample(nodes, self.sg.N)
            if self.savedyn:
                self.s_t.append(self.s.copy())
            dsNstep(smp)


    def run_c(self, steps: int | None = None, exName: str = ""):
        if steps is None:
            steps = self.simtime
        self.sg._export_edgel_bin(exName=exName)
        edgl = self.sg.path_exp_edgl
        cprog = LRGSG_CCORE_BIN / "contact_process"
        cmd = [str(cprog), f"{self.sg.N}", f"{self.mu}", f"{steps}", str(edgl)]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, check=True)
        self.s = np.frombuffer(result.stdout, dtype=np.int8)

    def run(self, **kwargs):
        if self.runlang.startswith("py"):
            self.run_py()
        elif self.runlang.startswith("C"):
            self.run_c(**kwargs)
        else:
            raise NotImplementedError("C implementation not available")