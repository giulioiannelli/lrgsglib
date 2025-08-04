from .common import *
from .BinDynSys import BinDynSys

class SignedRW(BinDynSys):
    dyn_UVclass = "signed_rw"
    id_string_signedrw = ""

    def __init__(self, sg: SignedGraph = Lattice2D, **kwargs) -> None:
        self.sg = sg
        self.dynpath = self.path_lrgsg
        super(BinDynSys, self).__init__(self.sg, **kwargs)

    def ds1step(self, nd: int):
        nodedict = dict(self.sg.H[nd])
        neighs = list(nodedict.keys())
        nn = np.random.choice(neighs)
        self.s[nd] = nodedict[nn]["weight"] * self.s[nn]

    def run_py(self):
        dsNstep = self.dsNstep()
        nodes = list(self.sg.H.nodes())
        for _ in range(self.simtime):
            smp = random.sample(nodes, self.sg.N)
            self.s_t.append(self.s.copy())
            # self.ds1step(np.random.randint(self.sg.N))
            dsNstep(smp)
    
    def run_c(self,
              adjfname: str = "",
              out_suffix: str = "",
              eqSTEP: int = 0):
        if eqSTEP:
            self.eqSTEP_def = eqSTEP
        if adjfname == "":
            adjfname = self.sg.std_fname
        out_suffix = out_suffix + self.id_string_voterdyn
        if out_suffix == "":
            out_suffix = '\'\''
        self.cprogram = [
            # pth_join(DIR_SRC, DIR_PCK, self.dyn_UVclass),
            # f"{DIR_SRC_DEFAULT}{DIR_PCK_DEFAULT}{self.dyn_UVclass}",
            f"{self.sg.N}",
            f"{self.sg.pflip}",
            f"{self.eqSTEP_def}",
            str(self.sg.path_data),
            adjfname,
            out_suffix
        ]
        subprocess.call(self.cprogram)

    def run(self, **kwargs):
        if self.runlang.startswith("py"):
            self.run_py()
        elif self.runlang.startswith("C"):
            self.run_c(**kwargs)