#include "LRGSG_rd.h"
#include "LRGSG_utils.h"
#include <string.h>
#include <math.h>

ReactionType parse_reaction_type(const char *s) {
    if (strsme(s, "fisher_kpp") || strsme(s, "fisher")) return RD_FISHER_KPP;
    if (strsme(s, "bistable"))  return RD_BISTABLE;
    if (strsme(s, "linear"))    return RD_LINEAR;
    if (strsme(s, "none"))      return RD_NONE;
    fprintf(stderr, "Unknown reaction type: %s\n", s);
    return RD_NONE;
}

static double reaction_term(double u, ReactionType type, double r, double a) {
    switch (type) {
    case RD_FISHER_KPP: return r * u * (1.0 - u);
    case RD_BISTABLE:   return u * (1.0 - u) * (u - a);
    case RD_LINEAR:     return r * u;
    case RD_NONE:       return 0.0;
    }
    return 0.0;
}

void rd_rhs(size_t N, const double *u, double *dudt, double t, void *params) {
    (void)t;
    RDParams *p = (RDParams *)params;

    for (size_t i = 0; i < N; ++i) {
        /* Diffusion: -D * L @ u = D * sum_j A_ij (u_j - u_i) */
        double diff = 0.0;
        size_t n_nn = p->neigh_len[i];
        for (size_t k = 0; k < n_nn; ++k) {
            size_t j = p->node_edges[i].neighbors[k];
            double w = p->node_edges[i].weights[k];
            diff += w * (u[j] - u[i]);
        }
        dudt[i] = p->D * diff + reaction_term(u[i], p->reaction, p->r, p->a);
    }
}
