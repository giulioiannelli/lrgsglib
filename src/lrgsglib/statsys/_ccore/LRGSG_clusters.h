#include "LRGSG_customs.h"

#ifndef __LRGSGCLUSTERSLIB_H_INC__
#define __LRGSGCLUSTERSLIB_H_INC__

/*
 * Incremental cluster-size-distribution tracker for single-spin spin dynamics
 * on a signed graph (shared _ccore; mirrors the Python reference
 * statsys/VoterModel/_cluster_tracker.py and is validated against its oracle).
 *
 * A *cluster* is a connected component of the subgraph of ACTIVE edges. The
 * activation predicate is selected by `rawspin`:
 *   rawspin = 0 ("satisfied") : edge (u,v) active iff sign(w_uv)*s_u*s_v > 0
 *   rawspin = 1 ("rawspin")   : edge (u,v) active iff s_u == s_v
 *
 * A single spin flip toggles every incident edge, so the partition is kept
 * incrementally: union-find-style MERGES with the now-agreeing neighbours are
 * O(1); a SPLIT can only happen when the node had >= 2 agreeing neighbours, and
 * only then runs a bounded local scan (cap `split_scan_cap`) that falls back to
 * an exact recompute of just that one affected component.
 *
 * The size distribution is recorded as a sparse histogram (size, count) pairs;
 * clusters_record() appends the current distribution to internal growable
 * storage, exposed flattened for the caller to marshal (pybind: list of dicts;
 * C subprocess: a file -- file wiring deferred, see the CP/cluster TODO).
 */

typedef struct ClusterCtx ClusterCtx;

/* Build the initial partition over the active-edge subgraph of `s`. The context
 * keeps a *reference* to `s`, `nlen` and `node_edges` (not copies), so the owner
 * must keep them alive and call clusters_flip() in lockstep with its flips. */
ClusterCtx *clusters_create(size_t N, const spin_tp s, size_tp nlen,
                            NodesEdges node_edges, int rawspin,
                            size_t split_scan_cap);

/* Update the partition for a flip of node `i`. MUST be called immediately BEFORE
 * the owner flips s[i] (the tracker reads the pre-flip configuration). */
void clusters_flip(ClusterCtx *ctx, size_t i);

/* Append the current cluster-size distribution (sparse histogram) as one record
 * to the internal growable storage. */
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
