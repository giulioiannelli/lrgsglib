#include "LRGSG_binsys.h"

#ifndef __LRGSGVMLIB_H_INC__
#define __LRGSGVMLIB_H_INC__

#define VTR_DIR "%s/voter/%s/"
#define SINI_FNAME VTR_DIR "s_" PSTR "%s" BINX
#define MAGN_FNAME VTR_DIR "m_" PSTR "%s" BINX
#define VTR_SOUT_FNAME VTR_DIR "sout_" PSTR "%s" BINX

/* Update rule (Axis A). Keep in sync with VoterModel/defaults.py RULE_CODE.
 * Signed substrate: the effective neighbour opinion is sign(w_ij) * s_j. */
typedef enum {
    VOTER_RULE_LINEAR    = 0,  /* copy a random neighbour                       */
    VOTER_RULE_MAJORITY  = 1,  /* sign(sum_j w_ij s_j), random tie-break        */
    VOTER_RULE_QVOTER    = 2,  /* q sampled w/ replacement; unanimous->copy     */
    VOTER_RULE_NONLINEAR = 3   /* P(+1)=f+^a/(f+^a+f-^a)                        */
} voter_rule_t;

/* Update schedule (Axis B). Keep in sync with VoterModel/defaults.py UPD_MODE_CODE. */
typedef enum {
    VOTER_UPD_ASYNC = 0,  /* N single-node updates, nodes drawn with replacement */
    VOTER_UPD_SYNC  = 1,  /* all nodes updated from a frozen snapshot            */
    VOTER_UPD_LINK  = 2   /* edge-update (rule='linear' only)                    */
} voter_upd_t;

/* Rule parameters (q-voter / nonlinear). */
typedef struct {
    voter_rule_t rule;
    size_t q;      /* q-voter: neighbours sampled with replacement */
    double eps;    /* q-voter: flip prob. when not unanimous       */
    double alpha;  /* nonlinear: exponent (alpha=1 -> linear)      */
} voter_params;

/* Apply the local rule to node `nd`, reading neighbour states from `s_read`;
 * returns the new spin (does not write). Serves both async (s_read=live) and
 * sync (s_read=snapshot) schedules. */
int8_t voter_apply_rule(size_t nd, const spin_tp s_read, size_tp nlen,
                        NodesEdges node_edges, voter_params p);

/* One asynchronous single-node update (in place). */
void voter_model_1step(size_t nd, spin_tp s, size_tp nlen,
                       NodesEdges node_edges, voter_params p);

/* One asynchronous sweep: N single-node updates, nodes drawn with replacement. */
void voter_model_Nstep(size_t N, spin_tp s, size_tp nlen,
                       NodesEdges node_edges, voter_params p);

/* One synchronous sweep: every node updated from snapshot `s` into `s_new`. */
void voter_sync_step(size_t N, spin_tp s, spin_tp s_new, size_tp nlen,
                     NodesEdges node_edges, voter_params p);

/* Build the cumulative-degree prefix sum cdeg[0..N]; returns total = 2*|E|. */
size_t voter_build_cdeg(size_t N, size_tp nlen, size_tp cdeg);

/* One link-update sweep: N edge-copies (linear-voter semantics). `cdeg`/`total`
 * come from voter_build_cdeg and give uniform directed-half-edge sampling. */
void voter_link_step(size_t N, spin_tp s, size_tp nlen, NodesEdges node_edges,
                     size_tp cdeg, size_t total);

/* Number of frustrated edges (w_ij s_i s_j < 0); 0 => absorbing configuration. */
size_t voter_count_frustrated(size_t N, spin_tp s, size_tp nlen,
                              NodesEdges node_edges);

#endif /* __LRGSGVMLIB_H_INC__ */
