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

/* Record the current configuration iff snap_out is requested AND either no
 * sample list is given (snap_sweeps == NULL => record every sweep) or the
 * current sweep `recorded` is the next sampled one. Sampled rows are packed
 * contiguously via `snap_row`, so a log-subsampled (or early-absorbing) run
 * allocates only the rows it actually keeps. Used identically by both impls. */
#define VOTER_SNAP_MAYBE()                                                     \
    do {                                                                       \
        if (snap_out && (snap_sweeps == NULL ||                               \
                (snap_si < n_snap_sweeps && snap_sweeps[snap_si] == recorded))) { \
            snap_record(&snap, &snap_cap_rows, snap_row, s, N);               \
            snap_row++;                                                        \
            if (snap_sweeps) snap_si++;                                       \
        }                                                                      \
    } while (0)

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
                                   const size_t *snap_sweeps, size_t n_snap_sweeps,
                                   size_t *snap_rows,
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
    size_t snap_row = 0, snap_si = 0;   /* packed snapshot row + sample cursor */
    while (recorded < n_sweeps) {
        if (total_f == 0) {                       /* frozen => absorbing */
            if (absorbing) {
                if (save_magn) magn[recorded] = calc_magn(N, s);
                VOTER_SNAP_MAYBE();
                if (track) clusters_record(cctx);
                *absorbed_at = (long)recorded;
                recorded++;
            } else {
                double m = save_magn ? calc_magn(N, s) : 0.0;
                while (recorded < n_sweeps) {      /* pad the frozen tail */
                    if (save_magn) magn[recorded] = m;
                    VOTER_SNAP_MAYBE();
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
            VOTER_SNAP_MAYBE();
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
    if (snap_rows) *snap_rows = snap_row;
    return recorded;
}

/* ================================================================== *
 * Composition-rejection (CR) sampler -- voter-specialised.
 *
 * Slepoy, Thompson & Plimpton, "A constant-time kinetic Monte Carlo
 * algorithm for simulation of large biochemical reaction networks",
 * J. Chem. Phys. 128, 205101 (2008): bucket events by rate so the
 * group is chosen by composition and the member by rejection, giving
 * O(1)-amortised event selection independent of N.
 *
 * Voter specialisation: node i flips at rate r_i = f_i / deg_i with the
 * integer frustration f_i in {0..deg_i}, so the rate spectrum is *exactly*
 * discrete. Bucketing on the integer key (deg_i, f_i) makes every node in a
 * bucket share one rate -> the within-bucket draw is UNIFORM with ZERO
 * rejection. Per event: choose a bucket prop. to count*rate by a linear scan
 * over occupied buckets (#buckets <= max_deg for a regular graph), draw a
 * member uniformly, flip, and re-bucket the node and its neighbours in O(1)
 * each -> O(deg) per event, dropping the Fenwick tree's log N factor.
 *
 * Buckets are kept in dynamic arrays with swap-remove (each node remembers its
 * position) and an occupied-key list (each key remembers its slot), so add /
 * remove / move are all O(1). total_rate is maintained incrementally,
 * eliminating the per-event prefix-sum. Recording / absorbing / snapshot /
 * cluster semantics are byte-for-byte the same as ctmc_run_impl.
 * ================================================================== */

typedef struct {
    size_t **gnode;    /* gnode[key] = dynamic array of node ids in bucket key */
    size_t  *gcap;     /* capacity of gnode[key]                               */
    size_t  *gcnt;     /* count in gnode[key]                                  */
    size_t  *npos;     /* per-node index within its bucket                     */
    size_t  *nkey;     /* per-node current key (0 = unbucketed, i.e. f_i == 0) */
    size_t  *occ;      /* list of currently-occupied keys                      */
    size_t  *occpos;   /* key -> slot in occ[] (n_keys = absent)               */
    size_t   n_occ;    /* number of occupied keys                              */
    size_t   n_keys;   /* (max_deg+1)^2; key = deg*stride + f                  */
    size_t   stride;   /* max_deg + 1                                          */
    double   total_rate;
} CRBuckets;

/* Rate f/deg recovered from the integer key (f = key % stride, deg = key/stride;
 * deg > 0 whenever a node is bucketed, since only f >= 1 nodes are stored). */
static inline double cr_key_rate(const CRBuckets *b, size_t key) {
    return (double)(key % b->stride) / (double)(key / b->stride);
}

/* Add node `nd` (degree `deg`, frustration `fi` >= 1) to its bucket. */
static void cr_add(CRBuckets *b, size_t nd, size_t deg, size_t fi) {
    size_t key = deg * b->stride + fi;
    if (b->gcnt[key] == b->gcap[key]) {
        size_t cap = b->gcap[key] ? b->gcap[key] << 1 : 8;
        size_t *p = realloc(b->gnode[key], cap * sizeof(*p));
        if (!p) { perror("realloc"); exit(EXIT_FAILURE); }
        b->gnode[key] = p;
        b->gcap[key] = cap;
    }
    if (b->gcnt[key] == 0) {                 /* bucket becomes occupied */
        b->occpos[key] = b->n_occ;
        b->occ[b->n_occ++] = key;
    }
    b->gnode[key][b->gcnt[key]] = nd;
    b->npos[nd] = b->gcnt[key];
    b->gcnt[key]++;
    b->nkey[nd] = key;
    b->total_rate += cr_key_rate(b, key);
}

/* Remove node `nd` from its bucket (no-op if currently unbucketed). */
static void cr_remove(CRBuckets *b, size_t nd) {
    size_t key = b->nkey[nd];
    if (key == 0) return;
    size_t pos = b->npos[nd];
    size_t last = b->gcnt[key] - 1;
    size_t moved = b->gnode[key][last];      /* swap the tail into the gap */
    b->gnode[key][pos] = moved;
    b->npos[moved] = pos;
    b->gcnt[key] = last;
    b->total_rate -= cr_key_rate(b, key);
    if (last == 0) {                         /* bucket emptied -> drop key */
        size_t oi = b->occpos[key];
        size_t mkey = b->occ[--b->n_occ];
        b->occ[oi] = mkey;
        b->occpos[mkey] = oi;
        b->occpos[key] = b->n_keys;          /* mark absent */
    }
    b->nkey[nd] = 0;
}

/* Move node `nd` (degree `deg`) to frustration `new_fi` (>= 0). */
static inline void cr_set(CRBuckets *b, size_t nd, size_t deg, size_t new_fi) {
    cr_remove(b, nd);
    if (new_fi > 0) cr_add(b, nd, deg, new_fi);
}

/* Draw a flipping node: bucket ~ count*rate (linear scan over occupied keys),
 * then a uniform member of that bucket (zero rejection). */
static size_t cr_sample_node(CRBuckets *b) {
    double target = RNG_dbl() * b->total_rate;
    size_t chosen = b->occ[b->n_occ - 1];    /* fp-rounding safe fallback */
    for (size_t o = 0; o < b->n_occ; o++) {
        size_t key = b->occ[o];
        double gsum = (double)b->gcnt[key] * cr_key_rate(b, key);
        if (target < gsum) { chosen = key; break; }
        target -= gsum;
    }
    size_t pick = (size_t)(RNG_u64() % b->gcnt[chosen]);
    return b->gnode[chosen][pick];
}

/*
 * CR core. Mirrors ctmc_run_impl exactly for the recording / absorbing / snap /
 * cluster bookkeeping (the `track` literal compile-folds the cluster branch the
 * same way); only the rate machinery (Fenwick -> (deg,f) buckets) and the node
 * draw differ. O(deg) per event, no log N.
 */
static inline size_t ctmc_run_cr_impl(size_t N, spin_tp s, size_tp nlen,
                                      NodesEdges node_edges, size_t n_sweeps,
                                      int save_magn, double *magn,
                                      spin_tp *snap_out, const size_t *snap_sweeps,
                                      size_t n_snap_sweeps, size_t *snap_rows,
                                      int absorbing, long *absorbed_at,
                                      const int track, ClusterCtx *cctx) {
    size_t max_deg = 0;
    for (size_t i = 0; i < N; i++) if (nlen[i] > max_deg) max_deg = nlen[i];
    size_t stride = max_deg + 1;
    size_t n_keys = stride * stride;         /* key = deg*stride + f < n_keys */

    size_tp f = __chMalloc(N * sizeof(*f));
    CRBuckets b;
    b.gnode  = calloc(n_keys, sizeof(*b.gnode));
    b.gcap   = calloc(n_keys, sizeof(*b.gcap));
    b.gcnt   = calloc(n_keys, sizeof(*b.gcnt));
    if (!b.gnode || !b.gcap || !b.gcnt) { perror("calloc"); exit(EXIT_FAILURE); }
    b.npos   = __chMalloc(N * sizeof(*b.npos));
    b.nkey   = __chMalloc(N * sizeof(*b.nkey));
    b.occ    = __chMalloc(n_keys * sizeof(*b.occ));
    b.occpos = __chMalloc(n_keys * sizeof(*b.occpos));
    for (size_t k = 0; k < n_keys; k++) b.occpos[k] = n_keys;  /* all absent */
    b.n_occ = 0; b.n_keys = n_keys; b.stride = stride; b.total_rate = 0.0;

    /* Initial frustration + buckets; total_f = 2 * frustrated edges. */
    long total_f = 0;
    for (size_t i = 0; i < N; i++) {
        size_t deg = nlen[i];
        size_t fi = deg ? node_frustration(i, s, deg, node_edges[i]) : 0;
        f[i] = fi;
        b.nkey[i] = 0;
        total_f += (long)fi;
        if (fi > 0) cr_add(&b, i, deg, fi);
    }

    double t = 0.0;          /* continuous time, in sweep units */
    size_t recorded = 0;     /* next integer sweep-time to record */
    spin_tp snap = NULL;     /* grown on demand, only when snap_out != NULL */
    size_t snap_cap_rows = 0;
    size_t snap_row = 0, snap_si = 0;   /* packed snapshot row + sample cursor */
    while (recorded < n_sweeps) {
        if (total_f == 0) {                       /* frozen => absorbing */
            if (absorbing) {
                if (save_magn) magn[recorded] = calc_magn(N, s);
                VOTER_SNAP_MAYBE();
                if (track) clusters_record(cctx);
                *absorbed_at = (long)recorded;
                recorded++;
            } else {
                double m = save_magn ? calc_magn(N, s) : 0.0;
                while (recorded < n_sweeps) {      /* pad the frozen tail */
                    if (save_magn) magn[recorded] = m;
                    VOTER_SNAP_MAYBE();
                    if (track) clusters_record(cctx);   /* distribution frozen too */
                    recorded++;
                }
            }
            break;
        }

        double dt = -log(1.0 - RNG_dbl()) / b.total_rate;  /* 1-u in (0,1] */
        double t_next = t + dt;

        /* State is constant on [t, t_next): emit any integer times in between. */
        while (recorded < n_sweeps && (double)recorded < t_next) {
            if (save_magn) magn[recorded] = calc_magn(N, s);
            VOTER_SNAP_MAYBE();
            if (track) clusters_record(cctx);
            recorded++;
        }
        if (recorded >= n_sweeps) break;

        /* Pick the flipping node ~ r_i / R (zero-rejection bucket draw). */
        size_t i = cr_sample_node(&b);
        size_t deg = nlen[i];
        s[i] = (int8_t)(-(int)s[i]);

        /* All edges at i toggle: new f_i = deg - old f_i. */
        size_t old_fi = f[i];
        size_t new_fi = deg - old_fi;
        f[i] = new_fi;
        total_f += (long)new_fi - (long)old_fi;
        cr_set(&b, i, deg, new_fi);

        /* Each neighbour j: edge (i,j) flipped status (uses the new s[i]). */
        NodeEdges ne = node_edges[i];
        for (size_t k = 0; k < deg; k++) {
            size_t j = ne.neighbors[k];
            int sign = (ne.weights[k] < 0.0) ? -1 : 1;
            int now_frustrated = (sign * (int)s[j] != (int)s[i]);
            int dj = now_frustrated ? 1 : -1;
            f[j] += (size_t)dj;            /* f[j] >= 1 before a -1 (shared edge) */
            total_f += dj;
            cr_set(&b, j, nlen[j], f[j]);
        }
        t = t_next;
    }

    for (size_t k = 0; k < n_keys; k++) free(b.gnode[k]);
    free(b.gnode); free(b.gcap); free(b.gcnt);
    free(b.npos); free(b.nkey); free(b.occ); free(b.occpos);
    free(f);
    if (snap_out) *snap_out = snap;
    if (snap_rows) *snap_rows = snap_row;
    return recorded;
}

/* Auto-select CR when the max degree is small enough that the (max_deg+1)^2
 * bucket table and the occupied-key scan stay cheap; above it the Fenwick BKL
 * (no degree-squared table, log-N selection) is the safer choice. */
#ifndef CTMC_CR_MAX_DEG
#define CTMC_CR_MAX_DEG 64
#endif

#undef VOTER_SNAP_MAYBE

size_t voter_ctmc_run(size_t N, spin_tp s, size_tp nlen, NodesEdges node_edges,
                      size_t n_sweeps, int save_magn, double *magn,
                      spin_tp *snap_out, const size_t *snap_sweeps,
                      size_t n_snap_sweeps, size_t *snap_rows,
                      int absorbing, long *absorbed_at, ClusterCtx *cctx) {
    *absorbed_at = -1;
    if (snap_out) *snap_out = NULL;
    if (snap_rows) *snap_rows = 0;
    if (n_sweeps == 0) return 0;

    /* Kernel choice: composition-rejection by default (O(deg)/event), with the
     * Fenwick BKL (O(deg log N)/event) as the dense-graph fallback. The
     * LRGSG_CTMC_KERNEL env var ("cr" | "fenwick" | "auto") forces one for
     * head-to-head benchmarking; unset/"auto" picks by max degree. */
    size_t max_deg = 0;
    for (size_t i = 0; i < N; i++) if (nlen[i] > max_deg) max_deg = nlen[i];
    int use_cr = (max_deg <= CTMC_CR_MAX_DEG);
    const char *force = getenv("LRGSG_CTMC_KERNEL");
    if (force) {
        if      (!strcmp(force, "cr"))      use_cr = 1;
        else if (!strcmp(force, "fenwick")) use_cr = 0;
    }

    /* Single runtime branch: pick the cluster-tracking or the plain loop ONCE.
     * Each *_impl instantiation has `track` as a literal, so its inner
     * `if (track)` checks compile away -- no per-iteration test either way. */
    if (use_cr) {
        if (cctx)
            return ctmc_run_cr_impl(N, s, nlen, node_edges, n_sweeps, save_magn,
                                    magn, snap_out, snap_sweeps, n_snap_sweeps, snap_rows, absorbing, absorbed_at, 1, cctx);
        return ctmc_run_cr_impl(N, s, nlen, node_edges, n_sweeps, save_magn,
                                magn, snap_out, snap_sweeps, n_snap_sweeps, snap_rows, absorbing, absorbed_at, 0, NULL);
    }
    if (cctx)
        return ctmc_run_impl(N, s, nlen, node_edges, n_sweeps, save_magn, magn,
                             snap_out, snap_sweeps, n_snap_sweeps, snap_rows, absorbing, absorbed_at, 1, cctx);
    return ctmc_run_impl(N, s, nlen, node_edges, n_sweeps, save_magn, magn,
                         snap_out, snap_sweeps, n_snap_sweeps, snap_rows, absorbing, absorbed_at, 0, NULL);
}
