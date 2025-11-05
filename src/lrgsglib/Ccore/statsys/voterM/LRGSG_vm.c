#include "LRGSG_vm.h"


/** 
 * @brief 
 * 
 * @param 
 * 
 * @return
 */
void voter_model_1step(size_t nd, spin_tp s, size_tp nlen, NodesEdges node_edges) {
    size_t degree = *(nlen + nd);
    if (!degree)
        return;
    size_t sel = (size_t)(RNG_u64() % degree);
    size_t neighbour = node_edges[nd].neighbors[sel];
    double weight = node_edges[nd].weights[sel];
    int sign = (weight < 0.0) ? -1 : 1;
    *(s + nd) = (int8_t)(sign * *(s + neighbour));
}

/** 
 * @brief 
 * 
 * @param 
 * 
 * @return
 */
void voter_model_Nstep(size_t N, spin_tp s, size_tp nlen, NodesEdges node_edges) {
    for (size_t i = 0; i < N; i++)
        voter_model_1step((size_t)(RNG_u64() % N), s, nlen, node_edges);
}