#include "LRGSG_vecdynsys.h"

#ifndef __LRGSG_POTTS_H_INC__
#define __LRGSG_POTTS_H_INC__

#define POTTS_DIR "%s/potts/%s/"
#define POTTS_SINI_FNAME POTTS_DIR "s_" PSTR "%s" BINX
#define POTTS_ENE_FNAME  POTTS_DIR "e_" PSTR "%s" BINX
#define POTTS_MAGN_FNAME POTTS_DIR "m_" PSTR "%s" BINX

/* Potts Metropolis: one sweep over all nodes */
void potts_metropolis_sweep(size_t N, int32_t *s, int q, double T,
                            size_t *neigh_len, NodeEdges *node_edges);

/* Calculate Potts energy: H = -sum A_ij delta(s_i, s_j) */
double calc_potts_energy(size_t N, const int32_t *s, size_t *neigh_len,
                         NodeEdges *node_edges);

#endif /* __LRGSG_POTTS_H_INC__ */
