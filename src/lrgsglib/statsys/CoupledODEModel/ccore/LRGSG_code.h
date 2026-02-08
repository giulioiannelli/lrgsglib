#include "LRGSG_contdynsys.h"

#ifndef __LRGSG_CODE_H_INC__
#define __LRGSG_CODE_H_INC__

#define CODE_DIR "%s/coupled_ode/%s/"
#define CODE_SINI_FNAME CODE_DIR "s_" PSTR "%s" BINX

/* Coupling type enum */
typedef enum {
    CODE_COUPLING_LINEAR  = 0,
    CODE_COUPLING_PRODUCT = 1,
} CouplingType;

/* Local function enum */
typedef enum {
    CODE_LOCAL_NONE     = 0,
    CODE_LOCAL_LINEAR   = 1,
    CODE_LOCAL_LOGISTIC = 2,
    CODE_LOCAL_FHN      = 3,
} LocalType;

/* Parameters struct */
typedef struct {
    size_t N;
    double coupling_strength;
    CouplingType coupling;
    LocalType local;
    double a;               /* linear coefficient */
    double r;               /* growth rate */
    double K;               /* carrying capacity */
    size_t *neigh_len;
    NodeEdges *node_edges;
} CODEParams;

void code_rhs(size_t N, const double *x, double *dxdt, double t, void *params);

CouplingType parse_coupling_type(const char *s);
LocalType parse_local_type(const char *s);

#endif /* __LRGSG_CODE_H_INC__ */
