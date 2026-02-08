#include "LRGSG_vecdynsys.h"

#ifndef __LRGSG_MULTISPEC_H_INC__
#define __LRGSG_MULTISPEC_H_INC__

#define MSPEC_DIR "%s/multi_species/%s/"
#define MSPEC_SINI_FNAME MSPEC_DIR "s_" PSTR "%s" BINX

/* Multi-species Metropolis: one sweep, flip one random component per node */
void multispec_metropolis_sweep(size_t N, int32_t *s, int k, int q, double T,
                                size_t *neigh_len, NodeEdges *node_edges);

#endif /* __LRGSG_MULTISPEC_H_INC__ */
