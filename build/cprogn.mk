# C filenames
FN_RBIMSIM0 = IsingSimulator0
FN_RBIMSIM1 = IsingSimulator1
FN_RBIMSIM2 = IsingSimulator2
FN_RBIMSIM3 = IsingSimulator3
FN_RBIMSIM4 = IsingSimulator4
FN_RBIMSIM5 = IsingSimulator5
FN_VMSIM0   = voter_model
FN_CPSIM0   = contact_process
FN_LRGSGLIB = LRGSG_utils sfmtrng 
SRC_BINDYNSYS = LRGSG_bindynsys
SRC_RBIM    = LRGSG_rbim
SRC_VM      = LRGSG_vm
SFMTSRC     = SFMT
#
FNS := $(FN_RBIMSIM0) $(FN_RBIMSIM1) $(FN_RBIMSIM2) \
       $(FN_RBIMSIM3) $(FN_RBIMSIM4) $(FN_RBIMSIM5) \
       $(FN_VMSIM0) $(FN_CPSIM0)
# only these go into PROGS
PROGS := $(addprefix $(LRGSG_CCORE_BIN)/, $(FNS))
# #
SRCCFILES.c     = $(addsuffix .c, $(FN_LRGSGLIB))
SRCCFILESBINDYNSYS.c = $(addsuffix .c, $(SRC_BINDYNSYS))
SRCCFILESRBIM.c = $(addsuffix .c, $(SRC_RBIM))
SRCCFILESVM.c   = $(addsuffix .c, $(SRC_VM))
SFMTFILES.c     = $(addsuffix .c, $(SFMTSRC))
# #
PATH_SRCC_FILES = $(addprefix $(LRGSG_LIB_CCORE)/, $(SRCCFILES.c))
PATH_SRCC_BINDYNSYS = $(addprefix $(LRGSG_CCORE_STATSYS)/, $(SRCCFILESBINDYNSYS.c))
PATH_SRCC_RBIM  = $(addprefix $(LRGSG_RBIM_SIMC)/, $(SRCCFILESRBIM.c))
PATH_SFMT_FILES = $(addprefix $(LRGSG_CCORE_SFMT)/, $(SFMTFILES.c))
PATH_SRCC_VM    = $(addprefix $(LRGSG_STATSYS_VM)/, $(SRCCFILESVM.c))
#