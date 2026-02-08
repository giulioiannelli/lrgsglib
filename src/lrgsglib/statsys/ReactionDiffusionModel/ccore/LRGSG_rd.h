#include "LRGSG_contdynsys.h"

#ifndef __LRGSG_RD_H_INC__
#define __LRGSG_RD_H_INC__

#define RD_DIR "%s/reaction_diffusion/%s/"
#define RD_SINI_FNAME RD_DIR "s_" PSTR "%s" BINX

/* Reaction type enum */
typedef enum {
    RD_FISHER_KPP = 0,
    RD_BISTABLE   = 1,
    RD_LINEAR     = 2,
    RD_NONE       = 3,
} ReactionType;

/* Parameters struct */
typedef struct {
    size_t N;
    double D;               /* diffusion constant */
    ReactionType reaction;
    double r;               /* growth rate (Fisher-KPP, linear) */
    double a;               /* threshold (bistable) */
    size_t *neigh_len;
    NodeEdges *node_edges;
} RDParams;

/* RHS function matching rhs_func_t signature */
void rd_rhs(size_t N, const double *u, double *dudt, double t, void *params);

/* Parse reaction type from string */
ReactionType parse_reaction_type(const char *s);

#endif /* __LRGSG_RD_H_INC__ */
