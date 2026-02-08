#include "LRGSG_kuramoto.h"
#include <math.h>

void kuramoto_rhs(size_t N, const double *theta, double *dtheta,
                  double t, void *params) {
    (void)t;
    KuramotoParams *p = (KuramotoParams *)params;

    for (size_t i = 0; i < N; ++i) {
        double coupling_sum = 0.0;
        size_t n_nn = p->neigh_len[i];
        for (size_t k = 0; k < n_nn; ++k) {
            size_t j = p->node_edges[i].neighbors[k];
            double w = p->node_edges[i].weights[k];
            coupling_sum += w * sin(theta[j] - theta[i]);
        }
        dtheta[i] = p->omega[i] + (p->K / (double)N) * coupling_sum;
    }
}

void kuramoto_step_rk4(size_t N, double *theta, double t, double dt,
                       KuramotoParams *par) {
    rk4_step(N, theta, kuramoto_rhs, t, dt, par);
    wrap_to_2pi(N, theta);
}
