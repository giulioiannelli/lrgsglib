#include "LRGSG_vecdynsys.h"

#ifndef __LRGSG_HEISENBERG_H_INC__
#define __LRGSG_HEISENBERG_H_INC__

#define HBERG_DIR "%s/heisenberg/%s/"
#define HBERG_SINI_FNAME HBERG_DIR "s_" PSTR "%s" BINX

/* Heisenberg Metropolis: one sweep over all nodes */
void heisenberg_metropolis_sweep(size_t N, double *spin, double T, double delta,
                                 size_t *neigh_len, NodeEdges *node_edges);

/* Heisenberg energy: H = -sum A_ij (n_i . n_j) */
double calc_heisenberg_energy(size_t N, const double *spin, size_t *neigh_len,
                              NodeEdges *node_edges);

/* Heisenberg magnetisation: |1/N sum n_i| */
double calc_heisenberg_magnetisation(size_t N, const double *spin);

#endif /* __LRGSG_HEISENBERG_H_INC__ */
