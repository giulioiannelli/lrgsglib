#include "LRGSG_vm.h"


/** 
 * @brief 
 * 
 * @param 
 * 
 * @return
 */
void voter_model_1step(size_t nd, spin_tp s, size_tp nlen, size_tp *neighs, double_p *edgl) {
    size_t sel_neigh_n = (size_t) (RNG_u64() % *(nlen + nd));
    size_t sel_neigh = *(*(neighs + nd) + sel_neigh_n);
    *(s + nd) = *(*(edgl + nd) + sel_neigh_n) * *(s + sel_neigh);
}

/** 
 * @brief 
 * 
 * @param 
 * 
 * @return
 */
void voter_model_Nstep(size_t N, spin_tp s, size_tp nlen, size_tp *neighs, double_p *edgl) {
    for (size_t i = 0; i < N; i++) 
        voter_model_1step((size_t) (RNG_u64() % N), s, nlen, neighs, edgl);
}