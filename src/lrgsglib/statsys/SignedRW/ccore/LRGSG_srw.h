#ifndef __LRGSG_SRW_H_INC__
#define __LRGSG_SRW_H_INC__

#include <stddef.h>
#include <stdint.h>

/* Stop reasons returned in stop_reason[w]. */
#define SRW_STOP_COVERED  0
#define SRW_STOP_DIED     1
#define SRW_STOP_TRAPPED  2

/* Absorbing-walker kernel (rule D3) with configurable X-node behaviour.
 *
 * Matches the Python reference in
 * lrgsglib.statsys.SignedRW._kernel._absorb_kernel.
 *
 * The signed adjacency is passed in CSR form together with two per-edge
 * boolean masks pre-computed on the Python side:
 *
 *   kill_mask[e]    == 1  ->  moving across edge e kills the walker
 *   reflect_mask[e] == 1  ->  moving across edge e leaves the walker in place
 *
 * Under 'reflect' x-node behaviour:
 *   kill_mask    = frust        (frustrated negative)
 *   reflect_mask = neg & ~frust (unfrustrated negative, i.e. X-node-incident)
 *
 * Under 'absorb' x-node behaviour:
 *   kill_mask    = neg          (any negative)
 *   reflect_mask = 0            (no reflections)
 *
 * The kernel advances every active walker by one step per outer iteration;
 * real-time `step` is the same as iteration count for the absorb rule.
 * It terminates when every walker is stopped (covered, died, or trapped).
 *
 * RNG: SFMT. The caller is responsible for seeding via `__set_seed_SFMT()`
 * (see statsys/_ccore/sfmtrng.h) before invoking this kernel.
 *
 * Output arrays must be pre-allocated by the caller. `visits_agg` is
 * zeroed by the kernel; the per-walker outputs are written directly.
 */
void run_absorb_walker(
    size_t N,
    size_t n_walkers,
    const size_t *neigh_ptr,            /* (N+1,) */
    const int64_t *neigh_indices,       /* (E,)   */
    const uint8_t *kill_mask,           /* (E,)   */
    const uint8_t *reflect_mask,        /* (E,)   */
    const int64_t *start_positions,     /* (n_walkers,) */
    int64_t stop_thresh,
    /* outputs */
    int64_t *visits_agg,                /* (N,)   */
    int64_t *unique_visits,             /* (n_walkers,) */
    int64_t *stop_step,                 /* (n_walkers,) */
    int8_t  *stop_reason,               /* (n_walkers,) */
    int64_t *final_position,            /* (n_walkers,) */
    size_t walkers_per_trial,           /* >0 = group contiguous walkers; 0 = disable */
    uint8_t *trial_bubble_out,          /* (n_trials*N,) or NULL; n_trials = n_walkers/walkers_per_trial */
    int track_unique                    /* 0 = skip per-walker bitmap + coverage check (fast path) */
);



/* Single-shot first-touch walker-count for two absorbing walkers on the same
 * signed CSR graph.  Alternates a walker from `node_a` and a walker from
 * `node_b`, growing two persistent bubbles.  Returns the smallest k >= 1 such
 * that after walker A's k-th run OR walker B's k-th run any visited site lies
 * in the other bubble.  If no intersection within `n_max` pairs, returns
 * `n_max + 1` (censored).  Each walker runs until death or `step_cap` steps.
 *
 * Semantics match `run_absorb_walker`'s per-step rule: reflect-aware advance,
 * record post-move visit, then resolve death.
 *
 * Outputs (may be NULL to skip):
 *   bubble_A_size_out / bubble_B_size_out = popcount of each bubble on return.
 *
 * Caller seeds SFMT via `__set_seed_SFMT()` (or the wrapper's `seed_rng`).
 */
int64_t run_absorb_walker_first_touch(
    size_t N,
    const size_t *neigh_ptr,
    const int64_t *neigh_indices,
    const uint8_t *kill_mask,
    const uint8_t *reflect_mask,
    int64_t node_a,
    int64_t node_b,
    size_t n_max,
    int64_t step_cap,
    int64_t *bubble_A_size_out,
    int64_t *bubble_B_size_out
);


#endif /* __LRGSG_SRW_H_INC__ */
