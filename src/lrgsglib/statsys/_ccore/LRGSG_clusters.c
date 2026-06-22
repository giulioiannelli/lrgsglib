#include "LRGSG_clusters.h"
#include <stdlib.h>
#include <stdio.h>

/* ------------------------------------------------------------------ *
 * Per-record full connected-components recompute over the active-edge
 * subgraph. No incremental partition state is kept: clusters_record()
 * runs one BFS sweep (O(N + E)) and emits the size histogram.
 *
 *   hist[size]    : scratch, number of components of each size (length N+1)
 *   vis/gen_vis   : generation-stamped BFS visited marks (O(1) clear)
 *   queue         : BFS frontier buffer
 *   rec_*         : recorded sparse histograms (record r is
 *                   rec_off[r]..rec_off[r+1] of rec_size[]/rec_cnt[])
 * ------------------------------------------------------------------ */
struct ClusterCtx {
    size_t N;
    const spin_tp s;          /* shared reference (owner mutates in lockstep) */
    size_tp nlen;
    NodesEdges ne;
    int rawspin;

    size_t *hist;             /* length N+1 */
    size_t *vis, gen_vis;
    size_t *queue;

    /* recorded sparse histograms: record r is rec_off[r]..rec_off[r+1] */
    size_t *rec_off, *rec_size, *rec_cnt;
    size_t n_rec, cap_off, n_ent, cap_ent;
};

static void *xmalloc(size_t n) {
    void *p = malloc(n);
    if (!p) { perror("malloc"); exit(EXIT_FAILURE); }
    return p;
}

/* active state of edge (u, ne[u].neighbors[k]) under the CURRENT s */
static inline int cl_active(const ClusterCtx *c, size_t u, size_t k) {
    size_t v = c->ne[u].neighbors[k];
    int su = (int)c->s[u], sv = (int)c->s[v];
    if (c->rawspin) return su == sv;
    int sign = (c->ne[u].weights[k] < 0.0) ? -1 : 1;
    return sign * su * sv > 0;
}

void clusters_record(ClusterCtx *c) {
    /* one connected-components pass; accumulate the size histogram in hist[] */
    size_t gv = ++c->gen_vis;
    size_t max_sz = 0;
    for (size_t start = 0; start < c->N; start++) {
        if (c->vis[start] == gv) continue;
        size_t sz = 0, qh = 0, qt = 0;
        c->vis[start] = gv;
        c->queue[qt++] = start;
        while (qh < qt) {
            size_t u = c->queue[qh++];
            sz++;
            size_t deg = c->nlen[u];
            for (size_t k = 0; k < deg; k++) {
                size_t v = c->ne[u].neighbors[k];
                if (c->vis[v] == gv) continue;
                if (!cl_active(c, u, k)) continue;
                c->vis[v] = gv;
                c->queue[qt++] = v;
            }
        }
        c->hist[sz]++;
        if (sz > max_sz) max_sz = sz;
    }

    /* emit the sparse histogram, clearing hist[] as we go (O(used sizes)) */
    if (c->n_rec + 1 >= c->cap_off) {
        c->cap_off = c->cap_off ? c->cap_off * 2 : 256;
        c->rec_off = realloc(c->rec_off, c->cap_off * sizeof(*c->rec_off));
        if (!c->rec_off) { perror("realloc"); exit(EXIT_FAILURE); }
    }
    for (size_t sz = 1; sz <= max_sz; sz++) {
        if (!c->hist[sz]) continue;
        if (c->n_ent >= c->cap_ent) {
            c->cap_ent = c->cap_ent ? c->cap_ent * 2 : 1024;
            c->rec_size = realloc(c->rec_size, c->cap_ent * sizeof(*c->rec_size));
            c->rec_cnt  = realloc(c->rec_cnt,  c->cap_ent * sizeof(*c->rec_cnt));
            if (!c->rec_size || !c->rec_cnt) { perror("realloc"); exit(EXIT_FAILURE); }
        }
        c->rec_size[c->n_ent] = sz;
        c->rec_cnt[c->n_ent]  = c->hist[sz];
        c->n_ent++;
        c->hist[sz] = 0;     /* reset for the next record */
    }
    c->n_rec++;
    c->rec_off[c->n_rec] = c->n_ent;
}

ClusterCtx *clusters_create(size_t N, const spin_tp s, size_tp nlen,
                            NodesEdges node_edges, int rawspin) {
    ClusterCtx *c = xmalloc(sizeof(*c));
    *(spin_tp *)&c->s = s;        /* init the const reference */
    c->N = N; c->nlen = nlen; c->ne = node_edges;
    c->rawspin = rawspin ? 1 : 0;

    c->hist  = calloc(N + 1, sizeof(size_t));
    c->vis   = calloc(N, sizeof(size_t));
    c->queue = xmalloc(N * sizeof(size_t));
    if (!c->hist || !c->vis) { perror("calloc"); exit(EXIT_FAILURE); }
    c->gen_vis = 0;

    c->rec_off = NULL; c->rec_size = NULL; c->rec_cnt = NULL;
    c->n_rec = 0; c->n_ent = 0; c->cap_ent = 0;
    c->cap_off = 256;
    c->rec_off = xmalloc(c->cap_off * sizeof(size_t));
    c->rec_off[0] = 0;       /* record 0 starts at entry 0 */
    return c;
}

size_t        clusters_num_records(const ClusterCtx *c) { return c->n_rec; }
const size_t *clusters_offsets(const ClusterCtx *c)     { return c->rec_off; }
const size_t *clusters_sizes(const ClusterCtx *c)       { return c->rec_size; }
const size_t *clusters_counts(const ClusterCtx *c)      { return c->rec_cnt; }

void clusters_free(ClusterCtx *c) {
    if (!c) return;
    free(c->hist); free(c->vis); free(c->queue);
    free(c->rec_off); free(c->rec_size); free(c->rec_cnt);
    free(c);
}
