# C filenames
FN_RBIMSIM0 = IsingSimulator0
FN_RBIMSIM1 = IsingSimulator1
FN_RBIMSIM2 = IsingSimulator2
FN_RBIMSIM3 = IsingSimulator3
FN_RBIMSIM4 = IsingSimulator4
FN_RBIMSIM5 = IsingSimulator5
FN_VMSIM0   = VoterSimulator0
FN_VMSIM1   = VoterSimulator1
FN_CPSIM0   = ContactSimulator0
FN_CPSIM1   = ContactSimulator1
FN_CPSIM1A  = ContactSimulator1a
FN_CPSIM1B  = ContactSimulator1b
FN_CPSIM1C  = ContactSimulator1c
FN_CPSIM1D  = ContactSimulator1d
FN_LRGSGLIB = LRGSG_utils sfmtrng
SRC_BINDYNSYS = LRGSG_bindynsys
SRC_RBIM    = LRGSG_rbim
SRC_VM      = LRGSG_vm
SRC_CP      = LRGSG_cp
SFMTSRC     = SFMT
#
FNS := $(FN_RBIMSIM0) $(FN_RBIMSIM1) $(FN_RBIMSIM2) \
       $(FN_RBIMSIM3) $(FN_RBIMSIM4) $(FN_RBIMSIM5) \
       $(FN_VMSIM0) $(FN_VMSIM1) $(FN_CPSIM0) $(FN_CPSIM1) \
       $(FN_CPSIM1A) $(FN_CPSIM1B) $(FN_CPSIM1C) $(FN_CPSIM1D)
# only these go into PROGS
PROGS := $(addprefix $(LRGSG_CCORE_BIN)/, $(FNS))
# #
SRCCFILES.c     = $(addsuffix .c, $(FN_LRGSGLIB))
SRCCFILESBINDYNSYS.c = $(addsuffix .c, $(SRC_BINDYNSYS))
SRCCFILESRBIM.c = $(addsuffix .c, $(SRC_RBIM))
SRCCFILESVM.c   = $(addsuffix .c, $(SRC_VM))
SRCCFILESCP.c   = $(addsuffix .c, $(SRC_CP))
SFMTFILES.c     = $(addsuffix .c, $(SFMTSRC))
# #
PATH_SRCC_FILES = $(addprefix $(LRGSG_LIB_CCORE)/, $(SRCCFILES.c))
PATH_SRCC_BINDYNSYS = $(addprefix $(LRGSG_CCORE_STATSYS)/, $(SRCCFILESBINDYNSYS.c))
PATH_SRCC_RBIM  = $(addprefix $(LRGSG_RBIM_SIMC)/, $(SRCCFILESRBIM.c))
PATH_SFMT_FILES = $(addprefix $(LRGSG_CCORE_SFMT)/, $(SFMTFILES.c))
PATH_SRCC_VM    = $(addprefix $(LRGSG_STATSYS_VM)/, $(SRCCFILESVM.c))
PATH_SRCC_CP    = $(addprefix $(LRGSG_STATSYS_CP)/, $(SRCCFILESCP.c))
#