#include "LRGSG_contdynsys.h"

#ifndef __LRGSG_KURAMOTO_H_INC__
#define __LRGSG_KURAMOTO_H_INC__

#define KUR_DIR "%s/kuramoto/%s/"
#define KUR_SINI_FNAME KUR_DIR "s_" PSTR "%s" BINX
#define KUR_HFLD_FNAME KUR_DIR "h_" PSTR "%s" BINX
#define KUR_ORDR_FNAME KUR_DIR "r_" PSTR "%s" BINX

/* Parameters struct passed through rhs_func_t void* */
typedef struct {
    size_t N;
    double K;               /* coupling strength */
    double *omega;           /* natural frequencies (length N) */
    size_t *neigh_len;
    NodeEdges *node_edges;
} KuramotoParams;

/* RHS function matching rhs_func_t signature */
void kuramoto_rhs(size_t N, const double *theta, double *dtheta,
                  double t, void *params);

/* One full Kuramoto integration step (RK4) */
void kuramoto_step_rk4(size_t N, double *theta, double t, double dt,
                       KuramotoParams *par);

#endif /* __LRGSG_KURAMOTO_H_INC__ */
