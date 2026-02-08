#include "LRGSG_binsys.h"

#ifndef __LRGSGVMLIB_H_INC__
#define __LRGSGVMLIB_H_INC__

#define VTR_DIR "%s/voter/%s/"
#define SINI_FNAME VTR_DIR "s_" PSTR "%s" BINX
#define MAGN_FNAME VTR_DIR "m_" PSTR "%s" BINX
#define VTR_SOUT_FNAME VTR_DIR "sout_" PSTR "%s" BINX

void voter_model_1step(size_t nd, spin_tp s, size_tp nlen, NodesEdges node_edges);
void voter_model_Nstep(size_t N, spin_tp s, size_tp nlen, NodesEdges node_edges);


#endif /* __LRGSGVMLIB_H_INC__ */