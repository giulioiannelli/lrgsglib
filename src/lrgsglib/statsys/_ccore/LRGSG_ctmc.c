#include "LRGSG_ctmc.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ *
 * Fenwick (binary-indexed) tree over node flip-rates, 1-indexed:
 * leaf position (i+1) holds rate[i]. point-update + prefix-sum +
 * find-by-cumulative all in O(log N), giving rate-proportional node
 * sampling for the rejection-free scheme.
 * ------------------------------------------------------------------ */
static inline size_t lowbit(size_t x) { return x & (~x + 1); }

/* Add `delta` to leaf for node `i` (0-indexed). */
static void fen_add(double *bit, size_t N, size_t i, double delta) {
    for (size_t x = i + 1; x <= N; x += lowbit(x))
        bit[x] += delta;
}

/* Prefix sum of leaves for nodes 0..i-1 (i.e. positions 1..i). */
static double fen_prefix(const double *bit, size_t i) {
    double s = 0.0;
    for (size_t x = i; x > 0; x -= lowbit(x)) s += bit[x];
    return s;
}

/* Smallest node k (0-indexed) with cumulative rate over 0..k strictly
 * exceeding `target`; assumes 0 <= target < total. Returns a node whose leaf
 * rate is > 0. */
static size_t fen_find(const double *bit, size_t N, double target) {
    size_t pos = 0, pw = 1;
    while ((pw << 1) <= N) pw <<= 1;
    for (; pw > 0; pw >>= 1) {
        if (pos + pw <= N && bit[pos + pw] <= target) {
            pos += pw;
            target -= bit[pos];
        }
    }
    return pos;   /* 0-indexed node = number of fully-accumulated leaves */
}

/* Count frustrated edges incident to node `nd` (sign(w)*s[nbr] != s[nd]). */
static size_t node_frustration(size_t nd, const spin_tp s, size_t deg,
                               NodeEdges ne) {
    size_t fi = 0;
    for (size_t k = 0; k < deg; k++) {
        int sign = (ne.weights[k] < 0.0) ? -1 : 1;
        if (sign * (int)s[ne.neighbors[k]] != (int)s[nd]) fi++;
    }
    return fi;
}

/* Append the current configuration as record `row` into a snapshot buffer that
 * grows on demand (capacity doubling) so an early-absorbing run allocates only
 * ~its record count, not the full n_sweeps*N. The caller frees *psnap. */
static inline void snap_record(spin_tp *psnap, size_t *pcap_rows, size_t row,
                               const spin_tp s, size_t N) {
    if (row + 1 > *pcap_rows) {
        size_t cap = *pcap_rows ? *pcap_rows : 64;
        while (cap < row + 1) cap <<= 1;
        spin_tp p = realloc(*psnap, cap * N * sizeof(*p));
        if (!p) { perror("realloc"); exit(EXIT_FAILURE); }
        *psnap = p;
        *pcap_rows = cap;
    }
    memcpy(*psnap + (size_t)row * N, s, N * sizeof(*s));
}

/*
 * Implementation core. `track` is a COMPILE-TIME literal (0 or 1) supplied by the
 * two call sites in voter_ctmc_run; being `static inline`, the compiler emits a
 * specialised copy per call site and constant-folds every `if (track)` away. So
 * the cluster branch is decided ONCE at dispatch -- there is no per-iteration
 * test in the no-tracking path (`cctx` is NULL there anyway).
 *
 * `snap_out` is an OPTIONAL full-trajectory sink: when non-NULL the kernel grows
 * a buffer on demand and copies the configuration at every recorded sweep time,
 * handing the buffer back via *snap_out (caller frees it). On-demand growth means
 * an early-absorbing run never allocates the full n_sweeps*N. Recording happens
 * only at the integer sweep times (O(records) memcpy), not per flip; NULL skips.
 */
static inline size_t ctmc_run_impl(size_t N, spin_tp s, size_tp nlen,
                                   NodesEdges node_edges, size_t n_sweeps,
                                   int save_magn, double *magn, spin_tp *snap_out,
                                   int absorbing, long *absorbed_at,
                                   const int track, ClusterCtx *cctx) {
    size_tp f    = __chMalloc(N * sizeof(*f));        /* frustrated edges/node */
    double *rate = __chMalloc(N * sizeof(*rate));     /* r_i = f_i / deg_i     */
    double *bit  = calloc(N + 1, sizeof(*bit));       /* Fenwick (zeroed)      */
    if (!bit) { perror("calloc"); exit(EXIT_FAILURE); }

    /* Initial frustration, rates and Fenwick; total_f = 2 * frustrated edges. */
    long total_f = 0;
    for (size_t i = 0; i < N; i++) {
        size_t deg = nlen[i];
        size_t fi = deg ? node_frustration(i, s, deg, node_edges[i]) : 0;
        f[i] = fi;
        total_f += (long)fi;
        double r = deg ? (double)fi / (double)deg : 0.0;
        rate[i] = r;
        if (r > 0.0) fen_add(bit, N, i, r);
    }

    double t = 0.0;          /* continuous time, in sweep units */
    size_t recorded = 0;     /* next integer sweep-time to record */
    spin_tp snap = NULL;     /* grown on demand, only when snap_out != NULL */
    size_t snap_cap_rows = 0;
    while (recorded < n_sweeps) {
        if (total_f == 0) {                       /* frozen => absorbing */
            if (absorbing) {
                if (save_magn) magn[recorded] = calc_magn(N, s);
                if (snap_out) snap_record(&snap, &snap_cap_rows, recorded, s, N);
                if (track) clusters_record(cctx);
                *absorbed_at = (long)recorded;
                recorded++;
            } else {
                double m = save_magn ? calc_magn(N, s) : 0.0;
                while (recorded < n_sweeps) {      /* pad the frozen tail */
                    if (save_magn) magn[recorded] = m;
                    if (snap_out) snap_record(&snap, &snap_cap_rows, recorded, s, N);
                    if (track) clusters_record(cctx);   /* distribution frozen too */
                    recorded++;
                }
            }
            break;
        }

        double R = fen_prefix(bit, N);
        double dt = -log(1.0 - RNG_dbl()) / R;     /* 1-u in (0,1] => finite */
        double t_next = t + dt;

        /* State is constant on [t, t_next): emit any integer times in between. */
        while (recorded < n_sweeps && (double)recorded < t_next) {
            if (save_magn) magn[recorded] = calc_magn(N, s);
            if (snap_out) snap_record(&snap, &snap_cap_rows, recorded, s, N);
            if (track) clusters_record(cctx);
            recorded++;
        }
        if (recorded >= n_sweeps) break;

        /* Pick the flipping node i ~ r_i / R and flip it. */
        size_t i = fen_find(bit, N, RNG_dbl() * R);
        if (i >= N) i = N - 1;
        size_t deg = nlen[i];
        s[i] = (int8_t)(-(int)s[i]);

        /* All edges at i toggle: new f_i = deg - old f_i. */
        size_t old_fi = f[i];
        size_t new_fi = deg - old_fi;
        f[i] = new_fi;
        total_f += (long)new_fi - (long)old_fi;
        double nr = deg ? (double)new_fi / (double)deg : 0.0;
        fen_add(bit, N, i, nr - rate[i]);
        rate[i] = nr;

        /* Each neighbour j: edge (i,j) flipped status (uses the new s[i]). */
        NodeEdges ne = node_edges[i];
        for (size_t k = 0; k < deg; k++) {
            size_t j = ne.neighbors[k];
            int sign = (ne.weights[k] < 0.0) ? -1 : 1;
            int now_frustrated = (sign * (int)s[j] != (int)s[i]);
            int dj = now_frustrated ? 1 : -1;
            f[j] += (size_t)dj;            /* f[j] >= 1 before a -1 (shared edge) */
            total_f += dj;
            double nrj = (double)f[j] / (double)nlen[j];   /* deg[j] >= 1 */
            fen_add(bit, N, j, nrj - rate[j]);
            rate[j] = nrj;
        }
        t = t_next;
    }

    free(f);
    free(rate);
    free(bit);
    if (snap_out) *snap_out = snap;
    return recorded;
}

size_t voter_ctmc_run(size_t N, spin_tp s, size_tp nlen, NodesEdges node_edges,
                      size_t n_sweeps, int save_magn, double *magn,
                      spin_tp *snap_out,
                      int absorbing, long *absorbed_at, ClusterCtx *cctx) {
    *absorbed_at = -1;
    if (snap_out) *snap_out = NULL;
    if (n_sweeps == 0) return 0;
    /* Single runtime branch: pick the cluster-tracking or the plain loop ONCE.
     * Each ctm_run_impl instantiation has `track` as a literal, so its inner
     * `if (track)` checks compile away -- no per-iteration test either way. */
    if (cctx)
        return ctmc_run_impl(N, s, nlen, node_edges, n_sweeps, save_magn, magn,
                             snap_out, absorbing, absorbed_at, 1, cctx);
    return ctmc_run_impl(N, s, nlen, node_edges, n_sweeps, save_magn, magn,
                         snap_out, absorbing, absorbed_at, 0, NULL);
}
