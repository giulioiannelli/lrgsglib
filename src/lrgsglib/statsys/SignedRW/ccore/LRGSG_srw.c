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
    int64_t *final_position,
    size_t walkers_per_trial,
    uint8_t *trial_bubble_out,
    int track_unique
) {
    /* per-walker bookkeeping */
    int64_t *pos = (int64_t *) malloc(n_walkers * sizeof(int64_t));
    uint8_t *alive = (uint8_t *) malloc(n_walkers * sizeof(uint8_t));
    uint8_t *stopped = (uint8_t *) malloc(n_walkers * sizeof(uint8_t));
    int64_t *step = (int64_t *) calloc(n_walkers, sizeof(int64_t));
    int64_t *unique = unique_visits;  /* alias: caller-owned buffer */

    /* Internal per-walker visited bitmap — required for per-walker
     * uniqueness tracking and the coverage-stop check.  Skipped entirely
     * on the fast path when track_unique == 0 (caller confirms no
     * coverage-based stopping is needed).  Never returned to caller. */
    uint8_t *visited = NULL;
    if (track_unique) {
        visited = (uint8_t *) calloc(n_walkers * N, sizeof(uint8_t));
    }

    /* Per-trial bubble buffer (OR of visited sites across walkers in the
     * same contiguous trial group).  The caller owns this buffer; we
     * zero-init it here if provided. */
    const bool want_bubbles =
        (trial_bubble_out != NULL) && (walkers_per_trial > 0);
    size_t n_trials = 0;
    if (want_bubbles) {
        if (n_walkers % walkers_per_trial != 0) {
            /* n_walkers must be an exact multiple of walkers_per_trial;
             * bail out without writing partial output. */
            free(pos); free(alive); free(stopped); free(step); free(visited);
            return;
        }
        n_trials = n_walkers / walkers_per_trial;
        memset(trial_bubble_out, 0, n_trials * N);
    }

    if (!pos || !alive || !stopped || !step || !unique ||
        (track_unique && !visited)) {
        free(pos); free(alive); free(stopped); free(step); free(visited);
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
        if (track_unique) {
            unique[w] = 1;
            visited[w * N + (size_t) p] = 1u;
        }
        visits_agg[p] += 1;
        if (want_bubbles) {
            size_t t = w / walkers_per_trial;
            trial_bubble_out[t * N + (size_t) p] = 1u;
        }
        stop_step[w] = -1;
        stop_reason[w] = -1;
        final_position[w] = p;
        if (track_unique && unique[w] >= stop_thresh) {
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
            if (track_unique) {
                uint8_t *bit = &visited[w * N + new_p];
                if (!*bit) {
                    *bit = 1u;
                    unique[w] += 1;
                }
            }
            visits_agg[new_p] += 1;
            if (want_bubbles) {
                size_t t = w / walkers_per_trial;
                trial_bubble_out[t * N + new_p] = 1u;
            }
            pos[w] = new_pos;
            final_position[w] = new_pos;

            /* resolve death (absorb-on-crossing) after the advance */
            if (is_kill) {
                alive[w] = 0u;
            }

            /* stop checks */
            if (track_unique && unique[w] >= stop_thresh) {
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
    /* unique and trial_bubble_out are caller-owned; do not free. */
    free(visited);
}


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
) {
    uint8_t *bubble_A = (uint8_t *) calloc(N, sizeof(uint8_t));
    uint8_t *bubble_B = (uint8_t *) calloc(N, sizeof(uint8_t));
    if (!bubble_A || !bubble_B) {
        free(bubble_A);
        free(bubble_B);
        return -1;
    }
    bubble_A[(size_t) node_a] = 1u;
    bubble_B[(size_t) node_b] = 1u;

    int64_t first_touch = -1;
    /* Trivial-start: if the start nodes coincide, touch is immediate (k=1). */
    if (node_a == node_b
        || bubble_A[(size_t) node_b]
        || bubble_B[(size_t) node_a]) {
        first_touch = 1;
    }

    for (size_t k = 1; k <= n_max && first_touch == -1; ++k) {
        /* walker A */
        int64_t pos = node_a;
        int alive = 1;
        for (int64_t s = 0; s < step_cap && alive && first_touch == -1; ++s) {
            size_t v = (size_t) pos;
            size_t beg = neigh_ptr[v];
            size_t end = neigh_ptr[v + 1];
            size_t deg = end - beg;
            if (deg == 0) break;
            size_t kk = (size_t) rand_u64_mod((uint64_t) deg);
            size_t e = beg + kk;
            int is_kill = (kill_mask[e] != 0);
            int is_refl = (reflect_mask[e] != 0);
            int64_t new_pos = pos;
            if (!is_refl) new_pos = neigh_indices[e];
            pos = new_pos;
            if (!bubble_A[(size_t) pos]) bubble_A[(size_t) pos] = 1u;
            if (bubble_B[(size_t) pos]) {
                first_touch = (int64_t) k;
                break;
            }
            if (is_kill) alive = 0;
        }
        if (first_touch != -1) break;

        /* walker B */
        pos = node_b;
        alive = 1;
        for (int64_t s = 0; s < step_cap && alive && first_touch == -1; ++s) {
            size_t v = (size_t) pos;
            size_t beg = neigh_ptr[v];
            size_t end = neigh_ptr[v + 1];
            size_t deg = end - beg;
            if (deg == 0) break;
            size_t kk = (size_t) rand_u64_mod((uint64_t) deg);
            size_t e = beg + kk;
            int is_kill = (kill_mask[e] != 0);
            int is_refl = (reflect_mask[e] != 0);
            int64_t new_pos = pos;
            if (!is_refl) new_pos = neigh_indices[e];
            pos = new_pos;
            if (!bubble_B[(size_t) pos]) bubble_B[(size_t) pos] = 1u;
            if (bubble_A[(size_t) pos]) {
                first_touch = (int64_t) k;
                break;
            }
            if (is_kill) alive = 0;
        }
    }

    if (bubble_A_size_out) {
        int64_t s = 0;
        for (size_t i = 0; i < N; ++i) s += bubble_A[i];
        *bubble_A_size_out = s;
    }
    if (bubble_B_size_out) {
        int64_t s = 0;
        for (size_t i = 0; i < N; ++i) s += bubble_B[i];
        *bubble_B_size_out = s;
    }
    free(bubble_A);
    free(bubble_B);
    return (first_touch == -1) ? (int64_t) (n_max + 1) : first_touch;
}

