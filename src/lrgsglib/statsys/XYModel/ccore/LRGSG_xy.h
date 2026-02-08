#include "LRGSG_contdynsys.h"
#include "LRGSG_vecdynsys.h"

#ifndef __LRGSG_XY_H_INC__
#define __LRGSG_XY_H_INC__

#define XY_DIR "%s/xy/%s/"
#define XY_SINI_FNAME XY_DIR "s_" PSTR "%s" BINX

/* XY Metropolis: one sweep over all nodes */
void xy_metropolis_sweep(size_t N, double *theta, double T, double delta,
                         size_t *neigh_len, NodeEdges *node_edges);

/* XY energy: H = -sum A_ij cos(theta_i - theta_j) */
double calc_xy_energy(size_t N, const double *theta, size_t *neigh_len,
                      NodeEdges *node_edges);

/* XY magnetisation: |1/N sum exp(i theta)| */
double calc_xy_magnetisation(size_t N, const double *theta);

#endif /* __LRGSG_XY_H_INC__ */
