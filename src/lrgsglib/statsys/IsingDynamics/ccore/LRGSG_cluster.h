/**
 * @file LRGSG_cluster.h
 * @brief Wolff and Swendsen-Wang cluster algorithms for Ising model.
 *
 * Supports signed (frustrated) graphs via coupling-aware bond activation:
 *   p_add = 1 - exp(-2*beta*|w|)  when s[i]*s[j]*sign(w) > 0
 *   p_add = 0                      otherwise
 */
#ifndef LRGSG_CLUSTER_H
#define LRGSG_CLUSTER_H

#include "LRGSG_binsys.h"

/**
 * Single Wolff cluster step.
 *
 * Picks a random seed, grows a cluster via BFS with bond probability
 * 1-exp(-2|w|/T), and flips the entire cluster.
 *
 * @param N     Number of spins.
 * @param T     Temperature (must be > 0).
 * @param s     Spin array (±1, modified in place).
 * @param nlen  Degree array.
 * @param ne    Adjacency structure.
 * @return      Size of the flipped cluster.
 */
size_t wolff_step(size_t N, double T, spin_tp s, size_tp nlen, NodesEdges ne);

/**
 * One Wolff "sweep": repeat wolff_step until ~N spins have been flipped.
 *
 * @return Total number of spins flipped in this sweep.
 */
size_t wolff_sweep(size_t N, double T, spin_tp s, size_tp nlen, NodesEdges ne);

/**
 * One Swendsen-Wang step.
 *
 * 1. Activate bonds between aligned spins with prob 1-exp(-2|w|/T).
 * 2. Find connected components via union-find.
 * 3. Flip each component with probability 0.5.
 *
 * @param N     Number of spins.
 * @param T     Temperature (must be > 0).
 * @param s     Spin array (±1, modified in place).
 * @param nlen  Degree array.
 * @param ne    Adjacency structure.
 * @return      Number of clusters found.
 */
size_t sw_step(size_t N, double T, spin_tp s, size_tp nlen, NodesEdges ne);

#endif /* LRGSG_CLUSTER_H */
