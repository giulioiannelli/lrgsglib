import networkx as nx
import numpy as np
import random_walker
import argparse

def generate_lattice(width, height):
    G = nx.grid_2d_graph(width, height)
    H = nx.convert_node_labels_to_integers(G)
    edge_list = list(H.edges())
    return edge_list, list(G.nodes())

def assign_random_signs(edge_list, p_flip=0.5):
    rand_p = np.random.rand(len(edge_list))
    signs = {edge: 1 if p_flip < rand_p[i] else -1 for i, edge in enumerate(edge_list)}
    return signs

def run_single_realization(edge_list, signs, width, height, coordination_number, pos):
    walker = random_walker.RandomWalker(edge_list, signs, width, height, coordination_number)
    start_node = width * (height // 2) + (width // 2)
    walker.walk(start_node)
    stopped_node = walker.get_stopped_node()
    original_stopped_node = pos[stopped_node]
    return walker.get_step_distance(), walker.get_euclidean_distance(), original_stopped_node

def run_multiple_realizations(width, height, pflip, navg, navg2, coordination_number):
    step_distances = []
    euclidean_distances = []
    stopped_nodes = []
    
    edge_list, pos = generate_lattice(width, height)
    for _ in range(navg2):
        signs = assign_random_signs(edge_list, pflip)
        
        avg_step_distance = 0
        avg_euclidean_distance = 0
        realization_stopped_nodes = []
        
        for _ in range(navg):
            step_distance, euclidean_distance, stopped_node = run_single_realization(edge_list, signs, width, height, coordination_number, pos)
            avg_step_distance += step_distance
            avg_euclidean_distance += euclidean_distance
            realization_stopped_nodes.append(stopped_node)
        
        avg_step_distance /= navg
        avg_euclidean_distance /= navg
        
        step_distances.append(avg_step_distance)
        euclidean_distances.append(avg_euclidean_distance)
        stopped_nodes.append(realization_stopped_nodes)
    
    avg_step_dist = sum(step_distances) / len(step_distances)
    avg_euclidean_dist = sum(euclidean_distances) / len(euclidean_distances)
    
    return avg_step_dist, avg_euclidean_dist, stopped_nodes

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Random Walker Simulation")
    parser.add_argument("--width", type=int, default=50, help="Width of the lattice")
    parser.add_argument("--height", type=int, default=50, help="Height of the lattice")
    parser.add_argument("--pflip", type=float, default=0.5, help="Probability of a negative sign")
    parser.add_argument("--navg", type=int, default=20, help="Number of realizations to average")
    parser.add_argument("--navg2", type=int, default=10, help="Number of times to repeat averaging")
    parser.add_argument("--coordination_number", type=int, default=4, help="Coordination number for the lattice")
    
    args = parser.parse_args()
    
    avg_step_dist, avg_euclidean_dist, stopped_nodes = run_multiple_realizations(
        args.width, args.height, args.pflip, args.navg, args.navg2, args.coordination_number)
    
    print("Average Step Distance:", avg_step_dist)
    print("Average Euclidean Distance:", avg_euclidean_dist)
    # print("Stopped Nodes:", stopped_nodes)
