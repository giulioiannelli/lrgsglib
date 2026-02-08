# C filenames
# ============================================================================
# Legacy Ising Simulators (old naming - kept for backwards compatibility)
FN_RBIMSIM0 = IsingSimulator0
FN_RBIMSIM1 = IsingSimulator1
FN_RBIMSIM2 = IsingSimulator2
FN_RBIMSIM3 = IsingSimulator3
FN_RBIMSIM4 = IsingSimulator4
FN_RBIMSIM5 = IsingSimulator5

# ============================================================================
# New Ising Simulator naming convention:
#   Number = Algorithm: 1=Metropolis, 2=Wolff(future), 3=SA, 4=PT, 5=SW(future)
#   Letter = Output: a=final, b=E/M, c=snapshots, d=clusters, e=eigvec, f=exchange
# ============================================================================

# Algorithm 1: Glauber-Metropolis
FN_ISINGSIM1B  = IsingSimulator1b

# Algorithm 3: Simulated Annealing
FN_ISINGSIM3B  = IsingSimulator3b

# Algorithm 4: Parallel Tempering
FN_ISINGSIM4B  = IsingSimulator4b

# Aggregate list of new Ising simulators
FNS_ISING_NEW := $(FN_ISINGSIM1B) $(FN_ISINGSIM3B) $(FN_ISINGSIM4B)

# Legacy Ising simulators aggregate
FNS_ISING_LEGACY := $(FN_RBIMSIM0) $(FN_RBIMSIM1) $(FN_RBIMSIM2) \
                    $(FN_RBIMSIM3) $(FN_RBIMSIM4) $(FN_RBIMSIM5)

# ============================================================================
# Voter Model Simulators
FN_VMSIM0   = VoterSimulator0
FN_VMSIM1   = VoterSimulator1

# Contact Process Simulators
# Unified simulator (consolidates 1, 1a-1g)
FN_CPSIM    = ContactSimulator
# Legacy variants (kept for backwards compatibility, deprecated)
FN_CPSIM0   = ContactSimulator0
FN_CPSIM1   = ContactSimulator1
FN_CPSIM1A  = ContactSimulator1a
FN_CPSIM1B  = ContactSimulator1b
FN_CPSIM1C  = ContactSimulator1c
FN_CPSIM1D  = ContactSimulator1d
FN_CPSIM1E  = ContactSimulator1e
FN_CPSIM1F  = ContactSimulator1f
FN_CPSIM1G  = ContactSimulator1g

# ============================================================================
# Continuous Dynamics Simulators
FN_KURSIM0    = KuramotoSimulator0
FN_RDSIM0     = ReactionDiffusionSimulator0
FN_CODESIM0   = CoupledODESimulator0

# Vector Dynamics Simulators
FN_POTTSSIM0  = PottsSimulator0
FN_XYSIM0     = XYSimulator0
FN_HBERGSIM0  = HeisenbergSimulator0
FN_MSPECSIM0  = MultiSpeciesSimulator0

# ============================================================================
# Shared Source Files
FN_LRGSGLIB = LRGSG_utils sfmtrng
SRC_BINDYNSYS    = LRGSG_bindynsys
SRC_CONTDYNSYS   = LRGSG_contdynsys
SRC_VECDYNSYS    = LRGSG_vecdynsys
SRC_RBIM    = LRGSG_rbim
SRC_SA      = LRGSG_sa
SRC_PT      = LRGSG_pt
SRC_VM      = LRGSG_vm
SRC_CP      = LRGSG_cp
SRC_KURAMOTO = LRGSG_kuramoto
SRC_RD       = LRGSG_rd
SRC_CODE     = LRGSG_code
SRC_POTTS    = LRGSG_potts
SRC_XY       = LRGSG_xy
SRC_HBERG    = LRGSG_heisenberg
SRC_MSPEC    = LRGSG_multispec
SFMTSRC     = SFMT

# ============================================================================
# All programs to build
# Note: FN_CPSIM (unified) replaces legacy FN_CPSIM1, FN_CPSIM1A-1G
# Legacy variants still built for backwards compatibility
FNS := $(FNS_ISING_LEGACY) $(FNS_ISING_NEW) \
       $(FN_VMSIM0) $(FN_VMSIM1) $(FN_CPSIM) $(FN_CPSIM0) $(FN_CPSIM1) \
       $(FN_CPSIM1A) $(FN_CPSIM1B) $(FN_CPSIM1C) $(FN_CPSIM1D) $(FN_CPSIM1E) $(FN_CPSIM1F) $(FN_CPSIM1G) \
       $(FN_KURSIM0) $(FN_RDSIM0) $(FN_CODESIM0) \
       $(FN_POTTSSIM0) $(FN_XYSIM0) $(FN_HBERGSIM0) $(FN_MSPECSIM0)
# Ising binaries go to per-model bin/ dirs; build PROGS list per-model
PROGS_ISING := $(addprefix $(LRGSG_STATSYS_ISING_BIN)/, $(FNS_ISING_LEGACY) $(FNS_ISING_NEW))
PROGS_VM    := $(addprefix $(LRGSG_STATSYS_VM_BIN)/, $(FN_VMSIM0) $(FN_VMSIM1))
PROGS_CP    := $(addprefix $(LRGSG_STATSYS_CP_BIN)/, $(FN_CPSIM) $(FN_CPSIM0) $(FN_CPSIM1) \
               $(FN_CPSIM1A) $(FN_CPSIM1B) $(FN_CPSIM1C) $(FN_CPSIM1D) $(FN_CPSIM1E) $(FN_CPSIM1F) $(FN_CPSIM1G))
PROGS_KUR   := $(addprefix $(LRGSG_STATSYS_KUR_BIN)/, $(FN_KURSIM0))
PROGS_RD    := $(addprefix $(LRGSG_STATSYS_RD_BIN)/, $(FN_RDSIM0))
PROGS_CODE  := $(addprefix $(LRGSG_STATSYS_CODE_BIN)/, $(FN_CODESIM0))
PROGS_POTTS := $(addprefix $(LRGSG_STATSYS_POTTS_BIN)/, $(FN_POTTSSIM0))
PROGS_XY    := $(addprefix $(LRGSG_STATSYS_XY_BIN)/, $(FN_XYSIM0))
PROGS_HBERG := $(addprefix $(LRGSG_STATSYS_HBERG_BIN)/, $(FN_HBERGSIM0))
PROGS_MSPEC := $(addprefix $(LRGSG_STATSYS_MSPEC_BIN)/, $(FN_MSPECSIM0))
PROGS := $(PROGS_ISING) $(PROGS_VM) $(PROGS_CP) \
         $(PROGS_KUR) $(PROGS_RD) $(PROGS_CODE) \
         $(PROGS_POTTS) $(PROGS_XY) $(PROGS_HBERG) $(PROGS_MSPEC)
# #
SRCCFILES.c     = $(addsuffix .c, $(FN_LRGSGLIB))
SRCCFILESBINDYNSYS.c = $(addsuffix .c, $(SRC_BINDYNSYS))
SRCCFILESRBIM.c = $(addsuffix .c, $(SRC_RBIM))
SRCCFILESSA.c   = $(addsuffix .c, $(SRC_SA))
SRCCFILESPT.c   = $(addsuffix .c, $(SRC_PT))
SRCCFILESVM.c   = $(addsuffix .c, $(SRC_VM))
SRCCFILESCP.c   = $(addsuffix .c, $(SRC_CP))
SFMTFILES.c     = $(addsuffix .c, $(SFMTSRC))
# #
PATH_SRCC_FILES = $(addprefix $(LRGSG_STATSYS_CCORE)/, $(SRCCFILES.c))
PATH_SRCC_BINDYNSYS = $(addprefix $(LRGSG_STATSYS_CCORE)/, $(SRCCFILESBINDYNSYS.c))
PATH_SRCC_RBIM  = $(addprefix $(LRGSG_STATSYS_ISING)/, $(SRCCFILESRBIM.c))
PATH_SRCC_SA    = $(addprefix $(LRGSG_STATSYS_ISING)/, $(SRCCFILESSA.c))
PATH_SRCC_PT    = $(addprefix $(LRGSG_STATSYS_ISING)/, $(SRCCFILESPT.c))
PATH_SFMT_FILES = $(addprefix $(LRGSG_CCORE_SFMT)/, $(SFMTFILES.c))
PATH_SRCC_VM    = $(addprefix $(LRGSG_STATSYS_VM)/, $(SRCCFILESVM.c))
PATH_SRCC_CP    = $(addprefix $(LRGSG_STATSYS_CP)/, $(SRCCFILESCP.c))
# Continuous / Vector dynamics shared sources
SRCCFILES_CONTDYNSYS.c = $(addsuffix .c, $(SRC_CONTDYNSYS))
SRCCFILES_VECDYNSYS.c  = $(addsuffix .c, $(SRC_VECDYNSYS))
PATH_SRCC_CONTDYNSYS = $(addprefix $(LRGSG_STATSYS_CCORE)/, $(SRCCFILES_CONTDYNSYS.c))
PATH_SRCC_VECDYNSYS  = $(addprefix $(LRGSG_STATSYS_CCORE)/, $(SRCCFILES_VECDYNSYS.c))
# Per-model shared sources
SRCCFILES_KUR.c   = $(addsuffix .c, $(SRC_KURAMOTO))
SRCCFILES_RD.c    = $(addsuffix .c, $(SRC_RD))
SRCCFILES_CODE.c  = $(addsuffix .c, $(SRC_CODE))
SRCCFILES_POTTS.c = $(addsuffix .c, $(SRC_POTTS))
SRCCFILES_XY.c    = $(addsuffix .c, $(SRC_XY))
SRCCFILES_HBERG.c = $(addsuffix .c, $(SRC_HBERG))
SRCCFILES_MSPEC.c = $(addsuffix .c, $(SRC_MSPEC))
PATH_SRCC_KUR   = $(addprefix $(LRGSG_STATSYS_KUR)/, $(SRCCFILES_KUR.c))
PATH_SRCC_RD    = $(addprefix $(LRGSG_STATSYS_RD)/, $(SRCCFILES_RD.c))
PATH_SRCC_CODE  = $(addprefix $(LRGSG_STATSYS_CODE)/, $(SRCCFILES_CODE.c))
PATH_SRCC_POTTS = $(addprefix $(LRGSG_STATSYS_POTTS)/, $(SRCCFILES_POTTS.c))
PATH_SRCC_XY    = $(addprefix $(LRGSG_STATSYS_XY)/, $(SRCCFILES_XY.c))
PATH_SRCC_HBERG = $(addprefix $(LRGSG_STATSYS_HBERG)/, $(SRCCFILES_HBERG.c))
PATH_SRCC_MSPEC = $(addprefix $(LRGSG_STATSYS_MSPEC)/, $(SRCCFILES_MSPEC.c))
#