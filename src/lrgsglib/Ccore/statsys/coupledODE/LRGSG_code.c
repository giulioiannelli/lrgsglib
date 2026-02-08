#include "LRGSG_code.h"
#include "LRGSG_utils.h"
#include <string.h>
#include <math.h>

CouplingType parse_coupling_type(const char *s) {
    if (strsme(s, "linear") || strsme(s, "diffusive")) return CODE_COUPLING_LINEAR;
    if (strsme(s, "product")) return CODE_COUPLING_PRODUCT;
    fprintf(stderr, "Unknown coupling type: %s, defaulting to linear\n", s);
    return CODE_COUPLING_LINEAR;
}

LocalType parse_local_type(const char *s) {
    if (strsme(s, "none"))      return CODE_LOCAL_NONE;
    if (strsme(s, "linear"))    return CODE_LOCAL_LINEAR;
    if (strsme(s, "logistic"))  return CODE_LOCAL_LOGISTIC;
    if (strsme(s, "fitzhugh") || strsme(s, "fhn")) return CODE_LOCAL_FHN;
    fprintf(stderr, "Unknown local type: %s, defaulting to none\n", s);
    return CODE_LOCAL_NONE;
}

static double local_func(double x, LocalType type, double a, double r, double K) {
    switch (type) {
    case CODE_LOCAL_NONE:     return 0.0;
    case CODE_LOCAL_LINEAR:   return a * x;
    case CODE_LOCAL_LOGISTIC: return r * x * (1.0 - x / K);
    case CODE_LOCAL_FHN:      return x - x * x * x / 3.0;
    }
    return 0.0;
}

void code_rhs(size_t N, const double *x, double *dxdt, double t, void *params) {
    (void)t;
    CODEParams *p = (CODEParams *)params;

    for (size_t i = 0; i < N; ++i) {
        double coupling_sum = 0.0;
        size_t n_nn = p->neigh_len[i];
        for (size_t k = 0; k < n_nn; ++k) {
            size_t j = p->node_edges[i].neighbors[k];
            double w = p->node_edges[i].weights[k];
            switch (p->coupling) {
            case CODE_COUPLING_LINEAR:
                coupling_sum += w * (x[j] - x[i]);
                break;
            case CODE_COUPLING_PRODUCT:
                coupling_sum += w * x[j] * x[i];
                break;
            }
        }
        dxdt[i] = local_func(x[i], p->local, p->a, p->r, p->K)
                 + p->coupling_strength * coupling_sum;
    }
}
