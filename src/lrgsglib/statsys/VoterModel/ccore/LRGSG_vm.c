#include "LRGSG_vm.h"
#include <math.h>

uint32_t *seed_rand;

/**
 * @brief Apply the local voter rule to node `nd`, reading from `s_read`.
 *
 * The signed substrate enters via the effective neighbour opinion
 * sigma_j = sign(w_ij) * s_read[j]. Serves both the asynchronous (s_read =
 * live array) and synchronous (s_read = frozen snapshot) schedules.
 */
int8_t voter_apply_rule(size_t nd, const spin_tp s_read, size_tp nlen,
                        NodesEdges node_edges, voter_params p) {
    size_t degree = *(nlen + nd);
    if (degree == 0) return *(s_read + nd);   /* isolated node: no change */
    NodeEdges ne = node_edges[nd];

    switch (p.rule) {
    case VOTER_RULE_LINEAR: {
        size_t sel = (size_t)(RNG_u64() % degree);
        int sign = (ne.weights[sel] < 0.0) ? -1 : 1;
        return (int8_t)(sign * (int)s_read[ne.neighbors[sel]]);
    }
    case VOTER_RULE_MAJORITY: {
        double h = 0.0;
        for (size_t i = 0; i < degree; i++)
            h += ne.weights[i] * (double)s_read[ne.neighbors[i]];
        if (h > 0.0) return (int8_t)1;
        if (h < 0.0) return (int8_t)(-1);
        return (int8_t)((RNG_u64() % 2) ? 1 : -1);   /* random tie-break */
    }
    case VOTER_RULE_QVOTER: {
        int8_t cur = s_read[nd];
        int first = 0, have = 0, unanimous = 1;
        for (size_t k = 0; k < p.q; k++) {
            size_t sel = (size_t)(RNG_u64() % degree);   /* with replacement */
            int sign = (ne.weights[sel] < 0.0) ? -1 : 1;
            int sigma = sign * (int)s_read[ne.neighbors[sel]];
            if (!have) { first = sigma; have = 1; }
            else if (sigma != first) { unanimous = 0; break; }
        }
        if (unanimous) return (int8_t)first;            /* copy the unanimous */
        if (RNG_dbl() < p.eps) return (int8_t)(-cur);   /* else flip w.p. eps */
        return cur;
    }
    case VOTER_RULE_NONLINEAR: {
        size_t nplus = 0;
        for (size_t i = 0; i < degree; i++) {
            int sign = (ne.weights[i] < 0.0) ? -1 : 1;
            if (sign * (int)s_read[ne.neighbors[i]] == 1) nplus++;
        }
        double fp = (double)nplus / (double)degree;
        double fm = 1.0 - fp;
        double num = pow(fp, p.alpha);
        double den = num + pow(fm, p.alpha);
        double pplus = (den > 0.0) ? (num / den) : 0.5;
        return (int8_t)((RNG_dbl() < pplus) ? 1 : -1);
    }
    }
    return *(s_read + nd);   /* unreachable: rule validated upstream */
}

/** Asynchronous in-place single-node update under `p`. */
void voter_model_1step(size_t nd, spin_tp s, size_tp nlen,
                       NodesEdges node_edges, voter_params p) {
    s[nd] = voter_apply_rule(nd, s, nlen, node_edges, p);
}

/** Asynchronous sweep: N single-node updates, nodes drawn with replacement. */
void voter_model_Nstep(size_t N, spin_tp s, size_tp nlen,
                       NodesEdges node_edges, voter_params p) {
    for (size_t i = 0; i < N; i++)
        voter_model_1step((size_t)(RNG_u64() % N), s, nlen, node_edges, p);
}

/** Synchronous sweep: every node updated from snapshot `s` into `s_new`. */
void voter_sync_step(size_t N, spin_tp s, spin_tp s_new, size_tp nlen,
                     NodesEdges node_edges, voter_params p) {
    for (size_t nd = 0; nd < N; nd++)
        s_new[nd] = voter_apply_rule(nd, s, nlen, node_edges, p);
}

/** Cumulative-degree prefix sum cdeg[0..N]; returns total = 2*|E|. */
size_t voter_build_cdeg(size_t N, size_tp nlen, size_tp cdeg) {
    cdeg[0] = 0;
    for (size_t i = 0; i < N; i++) cdeg[i + 1] = cdeg[i] + nlen[i];
    return cdeg[N];
}

/** Link-update sweep: N edge-copies (linear-voter semantics). */
void voter_link_step(size_t N, spin_tp s, size_tp nlen, NodesEdges node_edges,
                     size_tp cdeg, size_t total) {
    (void)nlen;
    if (total == 0) return;
    for (size_t step = 0; step < N; step++) {
        size_t k = (size_t)(RNG_u64() % total);   /* uniform directed half-edge */
        /* binary search: largest v with cdeg[v] <= k */
        size_t lo = 0, hi = N;
        while (lo + 1 < hi) {
            size_t mid = (lo + hi) / 2;
            if (cdeg[mid] <= k) lo = mid; else hi = mid;
        }
        size_t v = lo;
        size_t local = k - cdeg[v];
        size_t u = node_edges[v].neighbors[local];
        int sign = (node_edges[v].weights[local] < 0.0) ? -1 : 1;
        if (RNG_u64() % 2) s[v] = (int8_t)(sign * (int)s[u]);  /* v copies u */
        else               s[u] = (int8_t)(sign * (int)s[v]);  /* u copies v */
    }
}

/** Frustrated-edge count (w_ij s_i s_j < 0); 0 => absorbing configuration. */
size_t voter_count_frustrated(size_t N, spin_tp s, size_tp nlen,
                              NodesEdges node_edges) {
    size_t cnt = 0;
    for (size_t v = 0; v < N; v++) {
        NodeEdges ne = node_edges[v];
        size_t deg = nlen[v];
        for (size_t i = 0; i < deg; i++)
            if (ne.weights[i] * (double)s[v] * (double)s[ne.neighbors[i]] < 0.0)
                cnt++;
    }
    return cnt / 2;   /* each undirected edge counted twice */
}
