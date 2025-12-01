#ifndef __LRGSGCPLIB_H_INC__
#define __LRGSGCPLIB_H_INC__

#include "LRGSG_binsys.h"

#define CP_DIR "%s/cntct/%s/"
#define SINI_FNAME CP_DIR "s_" PSTR "_%s" BINX
#define SOUT_FNAME CP_DIR "sout_" PSTR "_%s" BINX
#define DENS_FNAME CP_DIR "dens_" PSTR "_%s" BINX

typedef struct cp_frontier_t cp_frontier_t;

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
void cp_lambda_init(double gamma, const spin_tp state, size_t N,
                    const NodesEdges node_edges, const size_tp neigh_len,
                    double *lambda);
void cp_lambda_update_neighbors(double gamma, size_t node, int delta_state,
                                const NodesEdges node_edges, const size_tp neigh_len,
                                double *lambda);
double cp_calculate_rate(int8_t state, double lambda, cp_activation_func_t activation_func);
void cp_rate_init(size_t N, const spin_tp state, const double *lambda,
                  cp_activation_func_t activation_func, double *rates, double *total_rate);
void cp_rate_update(size_t node, const spin_tp state, const double *lambda,
                    cp_activation_func_t activation_func, double *rates, double *total_rate);
void cp_update_rates_neighbors(size_t node, const NodesEdges node_edges, const size_tp neigh_len,
                               const spin_tp state, const double *lambda,
                               cp_activation_func_t activation_func, double *rates,
                               double *total_rate);
int cp_gillespie_select_event(size_t N, const double *rates, double total_rate, size_t *selected_node);
int cp_gillespie_select_event_frontier(const cp_frontier_t *frontier, const double *rates,
                                       double total_rate, size_t *selected_node);

typedef struct {
    double gamma;
    size_t N;
    spin_tp state;
    double *lambda;
    const NodesEdges node_edges;
    const size_tp neigh_len;
    cp_activation_func_t activation_func;
} cp_cached_sim_t;

size_t cp_run_sweep_cached(cp_cached_sim_t *sim, size_t *sum_ptr);

typedef struct {
    double gamma;
    size_t N;
    spin_tp state;
    double *lambda;
    double *rates;
    double total_rate;
    const NodesEdges node_edges;
    const size_tp neigh_len;
    cp_frontier_t *frontier;  // Optional frontier, used when selecting only active/boundary nodes
    cp_activation_func_t activation_func;
    int use_frontier_selector;
    int (*select_node)(const struct cp_gillespie_sim_t *sim, size_t *node_out);
} cp_gillespie_sim_t;

int cp_gillespie_pick_node(const cp_gillespie_sim_t *sim, size_t *node_out);
void cp_gillespie_set_selector(cp_gillespie_sim_t *sim);
void cp_gillespie_apply_flip(cp_gillespie_sim_t *sim, size_t node, int8_t old_state, int8_t new_state,
                             size_t *sum_ptr);

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

/* Frontier tracking for sparse simulations (three-list system) */
typedef enum {
    NODE_INACTIVE = 0,  // Not in any list (state=0, no active neighbors)
    NODE_BOUNDARY = 1,  // In boundary list (state=0, has ≥1 active neighbor)
    NODE_ACTIVE = 2     // In active list (state=1)
} cp_node_status_t;

typedef struct cp_frontier_t {
    size_t *active_list;            // List of active nodes (state=1)
    size_t *boundary_list;          // List of boundary nodes (state=0, has active neighbor)
    size_t *pos_active;             // Position lookup for active_list
    size_t *pos_boundary;           // Position lookup for boundary_list
    cp_node_status_t *node_status;  // Status of each node
    size_t *active_neighbor_count;  // Count of active neighbors for each node
    size_t N;                       // Number of nodes (for resets)
    size_t active_count;            // Number in active list
    size_t boundary_count;          // Number in boundary list
    int use_lambda_boundary;        // If true, only track boundary nodes with lambda>0
} cp_frontier_t;

/* Frontier management */
void cp_frontier_init(cp_frontier_t *frontier, size_t N);
void cp_frontier_free(cp_frontier_t *frontier);
void cp_frontier_build(cp_frontier_t *frontier, const spin_tp state, const double *lambda,
                       size_t N, const NodesEdges node_edges, const size_tp neigh_len);
void cp_frontier_node_activate(cp_frontier_t *frontier, size_t node,
                               const NodesEdges node_edges, const size_tp neigh_len,
                               const spin_tp state, const double *lambda);
void cp_frontier_node_deactivate(cp_frontier_t *frontier, size_t node,
                                 const NodesEdges node_edges, const size_tp neigh_len,
                                 const spin_tp state, const double *lambda);

typedef struct {
    double gamma;
    size_t N;
    spin_tp state;
    double *lambda;
    const NodesEdges node_edges;
    const size_tp neigh_len;
    cp_frontier_t *frontier;
    cp_activation_func_t activation_func;
} cp_frontier_sim_t;

typedef struct {
    size_t flips;
    size_t frontier_size;
} cp_frontier_sweep_result_t;

cp_frontier_sweep_result_t cp_run_frontier_sweep(cp_frontier_sim_t *sim, size_t *sum_ptr);
size_t cp_run_dense_frontier_sweep(cp_frontier_sim_t *sim, size_t *sum_ptr);

#endif /* __LRGSGCPLIB_H_INC__ */
