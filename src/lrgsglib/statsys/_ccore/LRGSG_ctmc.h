#include "LRGSG_binsys.h"
#include "LRGSG_clusters.h"

#ifndef __LRGSGCTMCLIB_H_INC__
#define __LRGSGCTMCLIB_H_INC__

/*
 * Rejection-free continuous-time MC (Gillespie / BKL) for the linear voter on a
 * signed graph. This is the shared CTMC driver; for now it is wired into the
 * VoterModel only (the Contact Process can adopt the same active-set/Fenwick
 * driver later -- see .agents/todo/2026-06-19_wire-ctmc-contact-process.md).
 *
 * Model. Each node attempts to copy a uniform signed neighbour at rate 1, so one
 * unit of continuous time equals one sweep (N attempts). A copy flips node i iff
 * the chosen edge is frustrated, hence node i flips as a Poisson process of rate
 *   r_i = f_i / deg_i,
 * where f_i is the number of frustrated edges incident to i
 * (sign(w_ij) * s_j != s_i). The next flip is sampled directly -- node i with
 * probability r_i / R, waiting time dt ~ Exp(R), R = sum_i r_i -- skipping the
 * null events of the rejection method. Absorption (R = 0, no frustrated edge)
 * is detected exactly from an integer frustration counter.
 *
 * Two interchangeable event-selection kernels (statistically identical, both
 * exact CTMC) are dispatched inside voter_ctmc_run:
 *   - composition-rejection (DEFAULT) -- buckets nodes by the integer key
 *     (deg_i, f_i); since f_i in {0..deg_i} the rate spectrum is exactly
 *     discrete, so the within-bucket draw is uniform with zero rejection.
 *     O(deg) per event, no log N. (Slepoy, Thompson & Plimpton, J. Chem. Phys.
 *     128, 205101 (2008), specialised to the voter's discrete rates.)
 *   - Fenwick (binary-indexed tree) BKL -- O(deg log N) per event; the
 *     fallback for dense graphs (large max degree).
 * voter_ctmc_run picks composition-rejection when the max degree is <=
 * CTMC_CR_MAX_DEG, else Fenwick; the env var LRGSG_CTMC_KERNEL
 * ("cr" | "fenwick" | "auto") overrides the choice for benchmarking.
 */

/*
 * Run the CTMC in place on `s` until the continuous time reaches `n_sweeps` or
 * the configuration freezes (absorbing). The magnetization is sampled at the
 * integer sweep times 0, 1, ..., (the state is piecewise-constant between flips)
 * into `magn` (caller-allocated, length >= n_sweeps) when `save_magn` is set.
 *
 * Returns t_run = number of sweep-times recorded (== n_sweeps unless it froze
 * early with `absorbing`). With `absorbing` set, freezing stops the run and
 * *absorbed_at receives the recorded sweep index of the frozen configuration;
 * otherwise *absorbed_at is left at -1 and the frozen tail is padded out to
 * n_sweeps. `magn` may be NULL when `save_magn` is 0.
 *
 * `snap_out` is an OPTIONAL full-trajectory sink: when non-NULL the kernel
 * allocates a buffer holding exactly the recorded configurations (row-major
 * (n_recorded, N) int8), grown on demand so an early-absorbing run never
 * allocates the full n_sweeps*N, and stores it in *snap_out; the CALLER owns it
 * and must free() it. The row count equals the return value. NULL skips it.
 *
 * `cctx` is an OPTIONAL cluster-size-distribution tracker: when non-NULL the
 * kernel records the current distribution at every integer sweep time; when NULL
 * the cluster machinery is skipped entirely (zero overhead -- only when
 * requested).
 */
size_t voter_ctmc_run(size_t N, spin_tp s, size_tp nlen, NodesEdges node_edges,
                      size_t n_sweeps, int save_magn, double *magn,
                      spin_tp *snap_out,
                      int absorbing, long *absorbed_at, ClusterCtx *cctx);

#endif /* __LRGSGCTMCLIB_H_INC__ */
