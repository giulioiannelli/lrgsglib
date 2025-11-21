#ifndef __LRGSGCPLIB_H_INC__
#define __LRGSGCPLIB_H_INC__

#include "LRGSG_binsys.h"

#define CP_DIR "%s/cntct/%s/"
#define SINI_FNAME CP_DIR "s_" PSTR "_%s" BINX
#define SOUT_FNAME CP_DIR "sout_" PSTR "_%s" BINX
#define DENS_FNAME CP_DIR "dens_" PSTR "_%s" BINX

/* Contact-process helper routines */
double cp_infection_rate(size_t node, spin_tp state, size_t degree, NodeEdges edges);
double cp_recovery_rate(size_t node, double mu, spin_tp state, size_t degree, NodeEdges edges);

typedef enum {
    CP_ACTIVATION_TANH = 0,
    CP_ACTIVATION_RELU = 1
} cp_activation_t;

/* Function pointer type for activation functions */
typedef double (*cp_activation_func_t)(double lambda);

cp_activation_t cp_activation_from_string(const char *name);
cp_activation_func_t cp_get_activation_function(cp_activation_t activation);
double cp_linear_input(double gamma, spin_tp state, size_t degree, NodeEdges edges);

/* Active border data structure for optimized contact process */
typedef struct {
    size_t *nodes;        // List of active border nodes
    size_t count;         // Current number of nodes in border
    size_t capacity;      // Allocated capacity
    int8_t *in_border;    // Fast lookup: is node i in border?
    size_t *position;     // Position of node i in the nodes array (valid only if in_border[i]==1)
} ActiveBorder;

/* Active border management functions */
void cp_init_active_border(ActiveBorder *ab, size_t N);
void cp_free_active_border(ActiveBorder *ab);
void cp_add_to_border(ActiveBorder *ab, size_t node);
void cp_remove_from_border(ActiveBorder *ab, size_t idx);
int cp_has_active_neighbor(const spin_tp state, size_t degree, const NodeEdges edges_node);
void cp_build_initial_border(ActiveBorder *ab, const spin_tp state, size_t N,
                             const NodesEdges node_edges, const size_tp neigh_len);
void cp_update_border_after_flip(ActiveBorder *ab, size_t node, int8_t new_state,
                                 const spin_tp state, size_t N,
                                 const NodesEdges node_edges, const size_tp neigh_len);

/* Simulation utility */
int cp_reached_absorbing_state(size_t sum, size_t N, size_t t, size_t steps);

#endif /* __LRGSGCPLIB_H_INC__ */
