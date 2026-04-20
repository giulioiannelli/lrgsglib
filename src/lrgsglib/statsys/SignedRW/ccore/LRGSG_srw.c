/* Absorbing-walker kernel (rule D3) — C port of Python _absorb_kernel.
 *
 * The Python reference lives in
 * lrgsglib.statsys.SignedRW._kernel._absorb_kernel.  See the header file
 * for the boolean-mask encoding of 'reflect' vs 'absorb' X-node modes.
 */

#include "LRGSG_srw.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "sfmtrng.h"


/* Uniform integer in [0, modulus).  For the small moduli used here
 * (degree z ~ 3-6) the modulo bias is negligible. */
static inline uint64_t rand_u64_mod(uint64_t modulus) {
    return SFMTrng_u64() % modulus;
}


void run_absorb_walker(
    size_t N,
    size_t n_walkers,
    const size_t *neigh_ptr,
    const int64_t *neigh_indices,
    const uint8_t *kill_mask,
    const uint8_t *reflect_mask,
    const int64_t *start_positions,
    int64_t stop_thresh,
    int64_t *visits_agg,
    int64_t *unique_visits,
    int64_t *stop_step,
    int8_t  *stop_reason,
    int64_t *final_position
) {
    /* per-walker bookkeeping */
    int64_t *pos = (int64_t *) malloc(n_walkers * sizeof(int64_t));
    uint8_t *alive = (uint8_t *) malloc(n_walkers * sizeof(uint8_t));
    uint8_t *stopped = (uint8_t *) malloc(n_walkers * sizeof(uint8_t));
    int64_t *step = (int64_t *) calloc(n_walkers, sizeof(int64_t));
    int64_t *unique = unique_visits;  /* alias: caller-owned buffer */

    /* per-walker visited bitmap: n_walkers * N bytes */
    uint8_t *visited = (uint8_t *) calloc(n_walkers * N, sizeof(uint8_t));

    if (!pos || !alive || !stopped || !step || !unique || !visited) {
        free(pos); free(alive); free(stopped); free(step);
        free(visited);
        return;
    }

    /* init aggregated visit counter */
    memset(visits_agg, 0, N * sizeof(int64_t));

    /* init per-walker state & initial position bookkeeping */
    for (size_t w = 0; w < n_walkers; ++w) {
        int64_t p = start_positions[w];
        pos[w] = p;
        alive[w] = 1u;
        stopped[w] = 0u;
        unique[w] = 1;
        visited[w * N + (size_t) p] = 1u;
        visits_agg[p] += 1;
        stop_step[w] = -1;
        stop_reason[w] = -1;
        final_position[w] = p;
        if (unique[w] >= stop_thresh) {
            stop_step[w] = 0;
            stop_reason[w] = SRW_STOP_COVERED;
            stopped[w] = 1u;
        }
    }

    /* trap check (mirror the Python guard): a walker whose starting
     * position has all outgoing edges reflective and none killing is
     * marked stopped with reason 2. */
    for (size_t w = 0; w < n_walkers; ++w) {
        if (stopped[w]) continue;
        size_t v = (size_t) pos[w];
        size_t beg = neigh_ptr[v];
        size_t end = neigh_ptr[v + 1];
        bool all_refl = (end > beg);
        bool any_kill = false;
        for (size_t e = beg; e < end; ++e) {
            if (!reflect_mask[e]) { all_refl = false; }
            if (kill_mask[e])     { any_kill = true; break; }
        }
        if (all_refl && !any_kill) {
            stop_step[w] = 0;
            stop_reason[w] = SRW_STOP_TRAPPED;
            stopped[w] = 1u;
        }
    }

    /* main loop: advance every active walker by one step until all are stopped */
    for (;;) {
        bool all_stopped = true;
        for (size_t w = 0; w < n_walkers; ++w) {
            if (stopped[w]) continue;
            all_stopped = false;

            size_t v = (size_t) pos[w];
            size_t beg = neigh_ptr[v];
            size_t end = neigh_ptr[v + 1];
            size_t deg = end - beg;

            /* Pick a uniformly random edge. */
            size_t k = (size_t) rand_u64_mod((uint64_t) deg);
            size_t e = beg + k;

            bool is_kill = (kill_mask[e] != 0);
            bool is_refl = (reflect_mask[e] != 0);

            step[w] += 1;

            /* Python semantics: advance FIRST (reflect-aware), record the
             * post-move visit, THEN resolve death.  A walker that crossed
             * a killing edge is recorded as having visited the target
             * before dying there. */
            int64_t new_pos = pos[w];
            if (!is_refl) {
                new_pos = neigh_indices[e];
            }

            size_t new_p = (size_t) new_pos;
            uint8_t *bit = &visited[w * N + new_p];
            if (!*bit) {
                *bit = 1u;
                unique[w] += 1;
            }
            visits_agg[new_p] += 1;
            pos[w] = new_pos;
            final_position[w] = new_pos;

            /* resolve death (absorb-on-crossing) after the advance */
            if (is_kill) {
                alive[w] = 0u;
            }

            /* stop checks */
            if (unique[w] >= stop_thresh) {
                stop_step[w] = step[w];
                stop_reason[w] = SRW_STOP_COVERED;
                stopped[w] = 1u;
            } else if (!alive[w]) {
                stop_step[w] = step[w];
                stop_reason[w] = SRW_STOP_DIED;
                stopped[w] = 1u;
            }
        }
        if (all_stopped) break;
    }

    free(pos);
    free(alive);
    free(stopped);
    free(step);
    /* unique is caller-owned; do not free */
    free(visited);
}
