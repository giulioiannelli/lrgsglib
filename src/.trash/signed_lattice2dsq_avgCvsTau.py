from lrgsglib.core import *
warnings.filterwarnings("ignore")
#
#
STEPS = 1500
SPIKE_THRESHOLD = 0.05
TAUSTOPEXP = 5
GEOMETRY = 'squared'
#
pow2_m, pow2_M, pow2_A = 10, 11, 16
ssize_range = range(pow2_m, pow2_M, 2)
asize_range = range(pow2_m, pow2_M, 1)[:len(ssize_range)]
#
lsN = np.array([2**i for i in ssize_range])
lsL = np.sqrt(lsN).astype('int')
lsA = np.array([1000 for i in range(len(ssize_range))])
#np.array([2**(pow2_A-i) for i in asize_range])
#
for iL, iN, iA in zip(lsL, lsN, lsA):
    #
    if not os.path.isdir((path := f"{datPath_l2d(GEOMETRY)}N={iN}_navg={iA}/")):
        os.makedirs(path)
    #
    lattice = Lattice2D(side1 = iL, geometry = GEOMETRY)
    lattice.lsp_selection([{'kind': 'lin', 'start': 0.089, 'stop': 0.096, 'num': 7, 'rsf':3}])#lattice.default_dict_lsp(num_at=7)
    print(lattice.lsp)
    for pflip in lattice.lsp:
        savename = lambda idstr : f"{path}p={pflip:{pflip_fmt}}_{idstr}{eBIN}"
        if os.path.exists(savename('Sm1')):
            continue
        Sm1File = open(savename('Sm1'), "wb")
        # CspeFile = open(savename('Cspe'), "wb")
        slspecFile = open(savename('slspec'), "wb")
        #
        for nr in tqdm(range(iA), desc=f"replicas for L={iL:3d}, p={pflip:6.3g}"):
            while True:
                SLRG_obj = SignedLaplacianAnalysis(#
                    sg = lattice,
                    pflip = pflip,
                    tauMex = TAUSTOPEXP,
                    steps = STEPS
                )
                SLRG_obj.flip_random_fract_edges()
                SLRG_obj.computeS()
                index = np.where(np.diff(SLRG_obj.Cspe) > SPIKE_THRESHOLD)
                if not index[0].size:
                    break
            #
            Sm1BytesArray = bytes(SLRG_obj.Sm1)
            # CspeBytesArray = bytes(SLRG_obj.Cspe)
            slspecBytesArray = bytes(SLRG_obj.slspectrum)
            #
            Sm1File.write(Sm1BytesArray)
            # CspeFile.write(CspeBytesArray)
            slspecFile.write(slspecBytesArray)
        Sm1File.close()
        # CspeFile.close()
        slspecFile.close()
