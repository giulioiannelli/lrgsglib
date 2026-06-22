#include "LRGSG_clusters.h"
#include <stdlib.h>
#include <stdio.h>

/* ------------------------------------------------------------------ *
 * Partition representation.
 *   comp_of[node]   : component id owning the node
 *   comp_size[id]   : number of nodes in component id
 *   comp_head[id]   : an arbitrary member of component id (ring entry point)
 *   nxt/prv[node]   : doubly-linked CIRCULAR ring per component (O(1) splice)
 *   hist[size]      : number of components of each size (the observable)
 *   free_ids/n_free : stack of unused component ids (<= N components ever)
 * Scratch (generation-stamped to avoid O(N) clears):
 *   vis/gen_vis     : BFS visited marks
 *   aux/gen_aux     : target marks (scan) / in-buffer marks (recompute)
 *   queue,membuf,fragbuf,Abuf,Bbuf : BFS / classification work buffers
 * ------------------------------------------------------------------ */
struct ClusterCtx {
    size_t N;
    const spin_tp s;          /* shared reference (owner mutates in lockstep) */
    size_tp nlen;
    NodesEdges ne;
    int rawspin;
    size_t split_cap;

    size_t *comp_of, *comp_size, *comp_head, *nxt, *prv;
    size_t *hist;             /* length N+1 */
    size_t *free_ids, n_free;

    size_t *vis, gen_vis;
    size_t *aux, gen_aux;
    size_t *queue, *membuf, *fragbuf, *Abuf, *Bbuf;

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

/* ------------------------------------------------------------------ *
 * Component primitives
 * ------------------------------------------------------------------ */
static size_t alloc_comp(ClusterCtx *c) { return c->free_ids[--c->n_free]; }
static void   free_comp(ClusterCtx *c, size_t id) { c->free_ids[c->n_free++] = id; }

/* Build a component from `nodes[0..len-1]` (len >= 1); updates hist. */
static size_t make_comp(ClusterCtx *c, const size_t *nodes, size_t len) {
    size_t id = alloc_comp(c);
    for (size_t k = 0; k < len; k++) {
        size_t u = nodes[k];
        c->comp_of[u] = id;
        c->nxt[u] = nodes[(k + 1) % len];
        c->prv[u] = nodes[(k + len - 1) % len];
    }
    c->comp_head[id] = nodes[0];
    c->comp_size[id] = len;
    c->hist[len]++;
    return id;
}

static void make_singleton(ClusterCtx *c, size_t i) {
    size_t id = alloc_comp(c);
    c->comp_of[i] = id;
    c->nxt[i] = i;
    c->prv[i] = i;
    c->comp_head[id] = i;
    c->comp_size[id] = 1;
    c->hist[1]++;
}

/* Unlink node i from its (size >= 2) component; remainder stays one component.
 * Does NOT touch comp_of[i] (caller re-homes i). */
static void remove_node(ClusterCtx *c, size_t i) {
    size_t id = c->comp_of[i], sz = c->comp_size[id];
    size_t pi = c->prv[i], ni = c->nxt[i];
    c->nxt[pi] = ni;
    c->prv[ni] = pi;
    if (c->comp_head[id] == i) c->comp_head[id] = ni;
    c->hist[sz]--;
    c->comp_size[id] = sz - 1;
    c->hist[sz - 1]++;
}

/* Merge component b into component a (a survives); updates hist + frees b. */
static void merge_into(ClusterCtx *c, size_t a, size_t b) {
    c->hist[c->comp_size[a]]--;
    c->hist[c->comp_size[b]]--;
    /* reassign b's members to a (walk b's ring) */
    size_t hb = c->comp_head[b], u = hb;
    do { c->comp_of[u] = a; u = c->nxt[u]; } while (u != hb);
    /* splice the two circular lists */
    size_t ha = c->comp_head[a];
    size_t at = c->prv[ha], bt = c->prv[hb];
    c->nxt[at] = hb; c->prv[hb] = at;
    c->nxt[bt] = ha; c->prv[ha] = bt;
    c->comp_size[a] += c->comp_size[b];
    free_comp(c, b);
    c->hist[c->comp_size[a]]++;
}

/* ------------------------------------------------------------------ *
 * Split handling
 * ------------------------------------------------------------------ */
/* Bounded BFS from A[0] over active edges excluding i. Returns 1 if all of A is
 * still mutually connected (no split); 0 if a split is detected or the cap is
 * hit (caller must recompute the component exactly). a_len >= 2. */
static int bounded_scan(ClusterCtx *c, size_t i, const size_t *A, size_t a_len) {
    size_t gv = ++c->gen_vis, ga = ++c->gen_aux;
    for (size_t t = 1; t < a_len; t++) c->aux[A[t]] = ga;   /* targets = A[1:] */
    size_t remaining = a_len - 1;
    size_t src = A[0];
    c->vis[src] = gv;
    size_t qh = 0, qt = 0, reached = 1;
    c->queue[qt++] = src;
    while (qh < qt) {
        size_t u = c->queue[qh++];
        size_t deg = c->nlen[u];
        for (size_t k = 0; k < deg; k++) {
            size_t v = c->ne[u].neighbors[k];
            if (v == i || c->vis[v] == gv) continue;
            if (!cl_active(c, u, k)) continue;
            c->vis[v] = gv;
            if (c->aux[v] == ga) { if (--remaining == 0) return 1; }
            if (++reached > c->split_cap) return 0;     /* cap -> recompute */
            c->queue[qt++] = v;
        }
    }
    return remaining == 0 ? 1 : 0;   /* exhausted: split iff targets remain */
}

/* Exact recompute: re-partition (old component minus i) into active-edge
 * fragments; i becomes a fresh singleton. */
static void recompute_split(ClusterCtx *c, size_t i, size_t id_old) {
    /* gather members of id_old except i */
    size_t sz = c->comp_size[id_old], m = 0, u = c->comp_head[id_old];
    for (size_t cnt = 0; cnt < sz; cnt++) {
        if (u != i) c->membuf[m++] = u;
        u = c->nxt[u];
    }
    c->hist[sz]--;
    free_comp(c, id_old);
    make_singleton(c, i);

    size_t ga = ++c->gen_aux, gv = ++c->gen_vis;
    for (size_t k = 0; k < m; k++) c->aux[c->membuf[k]] = ga;   /* in-buffer */
    for (size_t k = 0; k < m; k++) {
        size_t start = c->membuf[k];
        if (c->vis[start] == gv) continue;
        size_t fl = 0, qh = 0, qt = 0;
        c->vis[start] = gv;
        c->queue[qt++] = start;
        c->fragbuf[fl++] = start;
        while (qh < qt) {
            size_t w = c->queue[qh++];
            size_t deg = c->nlen[w];
            for (size_t kk = 0; kk < deg; kk++) {
                size_t v = c->ne[w].neighbors[kk];
                if (v == i || c->aux[v] != ga || c->vis[v] == gv) continue;
                if (!cl_active(c, w, kk)) continue;
                c->vis[v] = gv;
                c->queue[qt++] = v;
                c->fragbuf[fl++] = v;
            }
        }
        make_comp(c, c->fragbuf, fl);
    }
}

static void detach(ClusterCtx *c, size_t i, const size_t *A, size_t a_len) {
    size_t id = c->comp_of[i];
    if (c->comp_size[id] == 1) return;          /* i alone; re-homed by merge */
    if (a_len <= 1) {                            /* leaf/isolated: no split */
        remove_node(c, i);
        make_singleton(c, i);
        return;
    }
    if (bounded_scan(c, i, A, a_len)) {          /* connected: no split */
        remove_node(c, i);
        make_singleton(c, i);
        return;
    }
    recompute_split(c, i, id);                   /* genuine split / cap */
}

/* ------------------------------------------------------------------ *
 * Public API
 * ------------------------------------------------------------------ */
void clusters_flip(ClusterCtx *c, size_t i) {
    size_t deg = c->nlen[i];
    if (deg == 0) return;
    size_t a_len = 0, b_len = 0;
    for (size_t k = 0; k < deg; k++) {
        size_t v = c->ne[i].neighbors[k];
        if (cl_active(c, i, k)) c->Abuf[a_len++] = v;
        else                    c->Bbuf[b_len++] = v;
    }
    detach(c, i, c->Abuf, a_len);
    if (b_len) {                                 /* merge i with B-neighbours */
        size_t keep = c->comp_of[i];
        for (size_t k = 0; k < b_len; k++) {
            size_t cv = c->comp_of[c->Bbuf[k]];
            if (cv == keep) continue;
            if (c->comp_size[keep] >= c->comp_size[cv]) merge_into(c, keep, cv);
            else { merge_into(c, cv, keep); keep = cv; }
        }
    }
}

void clusters_record(ClusterCtx *c) {
    if (c->n_rec + 1 >= c->cap_off) {
        c->cap_off = c->cap_off ? c->cap_off * 2 : 256;
        c->rec_off = realloc(c->rec_off, c->cap_off * sizeof(*c->rec_off));
        if (!c->rec_off) { perror("realloc"); exit(EXIT_FAILURE); }
    }
    for (size_t sz = 1; sz <= c->N; sz++) {
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
    }
    c->n_rec++;
    c->rec_off[c->n_rec] = c->n_ent;
}

ClusterCtx *clusters_create(size_t N, const spin_tp s, size_tp nlen,
                            NodesEdges node_edges, int rawspin,
                            size_t split_scan_cap) {
    ClusterCtx *c = xmalloc(sizeof(*c));
    *(spin_tp *)&c->s = s;        /* init the const reference */
    c->N = N; c->nlen = nlen; c->ne = node_edges;
    c->rawspin = rawspin ? 1 : 0;
    c->split_cap = split_scan_cap ? split_scan_cap : 1;

    c->comp_of   = xmalloc(N * sizeof(size_t));
    c->comp_size = xmalloc(N * sizeof(size_t));
    c->comp_head = xmalloc(N * sizeof(size_t));
    c->nxt       = xmalloc(N * sizeof(size_t));
    c->prv       = xmalloc(N * sizeof(size_t));
    c->hist      = calloc(N + 1, sizeof(size_t));
    c->free_ids  = xmalloc(N * sizeof(size_t));
    c->vis       = calloc(N, sizeof(size_t));
    c->aux       = calloc(N, sizeof(size_t));
    c->queue     = xmalloc(N * sizeof(size_t));
    c->membuf    = xmalloc(N * sizeof(size_t));
    c->fragbuf   = xmalloc(N * sizeof(size_t));
    c->Abuf      = xmalloc(N * sizeof(size_t));
    c->Bbuf      = xmalloc(N * sizeof(size_t));
    if (!c->hist || !c->vis || !c->aux) { perror("calloc"); exit(EXIT_FAILURE); }
    c->gen_vis = 0; c->gen_aux = 0;

    /* all component ids start free */
    for (size_t id = 0; id < N; id++) c->free_ids[id] = id;
    c->n_free = N;

    /* initial partition: connected components over the active-edge subgraph */
    size_t gv = ++c->gen_vis;
    for (size_t start = 0; start < N; start++) {
        if (c->vis[start] == gv) continue;
        size_t fl = 0, qh = 0, qt = 0;
        c->vis[start] = gv;
        c->queue[qt++] = start;
        c->fragbuf[fl++] = start;
        while (qh < qt) {
            size_t u = c->queue[qh++];
            size_t deg = nlen[u];
            for (size_t k = 0; k < deg; k++) {
                size_t v = node_edges[u].neighbors[k];
                if (c->vis[v] == gv) continue;
                if (!cl_active(c, u, k)) continue;
                c->vis[v] = gv;
                c->queue[qt++] = v;
                c->fragbuf[fl++] = v;
            }
        }
        make_comp(c, c->fragbuf, fl);
    }

    c->rec_off = NULL; c->rec_size = NULL; c->rec_cnt = NULL;
    c->n_rec = 0; c->cap_off = 0; c->n_ent = 0; c->cap_ent = 0;
    /* seed offsets with rec_off[0] = 0 */
    c->cap_off = 256;
    c->rec_off = xmalloc(c->cap_off * sizeof(size_t));
    c->rec_off[0] = 0;
    return c;
}

size_t        clusters_num_records(const ClusterCtx *c) { return c->n_rec; }
const size_t *clusters_offsets(const ClusterCtx *c)     { return c->rec_off; }
const size_t *clusters_sizes(const ClusterCtx *c)       { return c->rec_size; }
const size_t *clusters_counts(const ClusterCtx *c)      { return c->rec_cnt; }

void clusters_free(ClusterCtx *c) {
    if (!c) return;
    free(c->comp_of); free(c->comp_size); free(c->comp_head);
    free(c->nxt); free(c->prv); free(c->hist); free(c->free_ids);
    free(c->vis); free(c->aux); free(c->queue); free(c->membuf);
    free(c->fragbuf); free(c->Abuf); free(c->Bbuf);
    free(c->rec_off); free(c->rec_size); free(c->rec_cnt);
    free(c);
}
