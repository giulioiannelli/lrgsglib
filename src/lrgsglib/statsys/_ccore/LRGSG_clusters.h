#include "LRGSG_customs.h"

#ifndef __LRGSGCLUSTERSLIB_H_INC__
#define __LRGSGCLUSTERSLIB_H_INC__

/*
 * Cluster-size-distribution observable for spin dynamics on a signed graph
 * (shared _ccore; mirrors the Python reference
 * utils/statsys.py and is validated against it).
 *
 * A *cluster* is a connected component of the subgraph of ACTIVE edges. The
 * activation predicate is selected by `rawspin`:
 *   rawspin = 0 ("satisfied") : edge (u,v) active iff sign(w_uv)*s_u*s_v > 0
 *   rawspin = 1 ("rawspin")   : edge (u,v) active iff s_u == s_v
 *
 * clusters_record() recomputes the whole partition from scratch -- a single
 * connected-components pass over the active-edge subgraph (O(N + E)) -- and
 * appends the resulting size distribution. The owner calls it once per recorded
 * sweep; with ~N single-spin events per sweep this matches the cost of
 * incremental maintenance while remaining trivially correct (no per-flip hook is
 * needed, so the kernel's inner loop carries zero cluster overhead).
 *
 * The distribution is stored as a sparse histogram of (size, count) pairs,
 * exposed flattened for the caller to marshal (pybind: list of dicts;
 * C subprocess: a file -- file wiring deferred, see the CP/cluster TODO).
 */

typedef struct ClusterCtx ClusterCtx;

/* Allocate the recompute context. It keeps a *reference* to `s`, `nlen` and
 * `node_edges` (not copies), so the owner must keep them alive for its lifetime.
 * `rawspin` selects the activation predicate (see above). */
ClusterCtx *clusters_create(size_t N, const spin_tp s, size_tp nlen,
                            NodesEdges node_edges, int rawspin);

/* Re-point the context at state buffer `s`, for owners that double-buffer and
 * swap `s` between sweeps (e.g. the synchronous voter). The next
 * clusters_record() reads from this buffer. */
void clusters_set_state(ClusterCtx *ctx, const spin_tp s);

/* Recompute the active-edge connected components under the CURRENT `s` and
 * append the resulting cluster-size distribution as one record. */
void clusters_record(ClusterCtx *ctx);

/* Number of recorded distributions so far. */
size_t clusters_num_records(const ClusterCtx *ctx);

/* Flattened access to the recorded distributions. Record r occupies the half-open
 * range [offsets()[r], offsets()[r+1]) of sizes()/counts(); i.e. for that record
 * there are counts()[k] clusters of size sizes()[k]. `offsets` has length
 * num_records()+1. Pointers are owned by the context (valid until free). */
const size_t *clusters_offsets(const ClusterCtx *ctx);
const size_t *clusters_sizes(const ClusterCtx *ctx);
const size_t *clusters_counts(const ClusterCtx *ctx);

void clusters_free(ClusterCtx *ctx);

#endif /* __LRGSGCLUSTERSLIB_H_INC__ */
