#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

typedef struct {
    uint64_t u, v;
    double w;
} Edge;

// Function to read edges from a binary file
Edge* read_edges(const char *filename, int *edge_count) {
    FILE *file = fopen(filename, "rb");
    if (!file) {
        perror("Cannot open file");
        exit(EXIT_FAILURE);
    }

    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    fseek(file, 0, SEEK_SET);

    *edge_count = file_size / sizeof(Edge);
    Edge *edges = (Edge*)malloc(file_size);
    if (!edges) {
        perror("Memory allocation failed");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    fread(edges, sizeof(Edge), *edge_count, file);
    fclose(file);

    return edges;
}

// Function to create bidirectional edges and initialize neighs and neigh_len arrays
void process_edges(const char *filename, int N, size_t ***neighs, size_t **neigh_len) {
    int edge_count, new_edge_count;
    printf("Reading edges from %s\n", filename);
    Edge *edges = read_edges(filename, &edge_count);
    printf("Edge count: %d\n", edge_count);
    new_edge_count = 2 * edge_count;

    Edge *all_edges = (Edge*)malloc(new_edge_count * sizeof(Edge));
    if (!all_edges) {
        perror("Memory allocation failed");
        free(edges);
        exit(EXIT_FAILURE);
    }
    printf("New edge count: %d\n", new_edge_count);

    // Allocate memory for neigh_len and initialize to 0
    *neigh_len = (size_t*)calloc(N, sizeof(size_t));
    if (!*neigh_len) {
        perror("Memory allocation failed");
        free(edges);
        free(all_edges);
        exit(EXIT_FAILURE);
    }
    printf("Allocated memory for neigh_len\n");


    for (int i = 0; i < edge_count; i++) {
    // Check for out-of-bounds access
        if (edges[i].u >= N || edges[i].v >= N) {
            fprintf(stderr, "Error: node index out of bounds. u: %ld, v: %ld, N: %d\n", edges[i].u, edges[i].v, N);
            free(edges);
            free(all_edges);
            free(*neigh_len);
            exit(EXIT_FAILURE);
        }

        *(all_edges + 2 * i) = *(edges + i);
        all_edges[2 * i + 1].u = edges[i].v;
        all_edges[2 * i + 1].v = edges[i].u;
        all_edges[2 * i + 1].w = edges[i].w;

        (*neigh_len)[edges[i].u]++;
        (*neigh_len)[edges[i].v]++;
        
        printf("Edge %d: %ld %ld %f\n", i, edges[i].u, edges[i].v, edges[i].w);
    }

    // Allocate memory for neighs arrays based on neigh_len
    *neighs = (size_t**)malloc(N * sizeof(size_t*));
    if (!*neighs) {
        perror("Memory allocation failed");
        free(edges);
        free(all_edges);
        free(*neigh_len);
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < N; i++) {
        (*neighs)[i] = (size_t*)malloc((*neigh_len)[i] * sizeof(size_t));
        if (!(*neighs)[i]) {
            perror("Memory allocation failed");
            for (int j = 0; j < i; j++) {
                free((*neighs)[j]);
            }
            free(*neighs);
            free(edges);
            free(all_edges);
            free(*neigh_len);
            exit(EXIT_FAILURE);
        }
        (*neigh_len)[i] = 0; // Reset to use as an index
    }

    // Fill neighs arrays
    for (int i = 0; i < new_edge_count; i++) {
        int u = all_edges[i].u;
        int v = all_edges[i].v;
        (*neighs)[u][(*neigh_len)[u]++] = v;
    }

    free(edges);
    free(all_edges);
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Usage: %s <number_of_nodes> <edge_file>\n", argv[0]);
        return EXIT_FAILURE;
    }
    
    int N = atoi(argv[1]);
    printf("N: %d, dir: %s\n", N, argv[2]);
    size_t **neighs;
    size_t *neigh_len;

    process_edges(argv[2], N, &neighs, &neigh_len);

    // Example: Print neighbors for each node
    for (int i = 0; i < N; i++) {
        printf("Neighbors of node %d: ", i);
        for (int j = 0; j < neigh_len[i]; j++) {
            printf("%zu ", neighs[i][j]);
        }
        printf("\n");
    }

    // Free memory
    for (int i = 0; i < N; i++) {
        free(neighs[i]);
    }
    free(neighs);
    free(neigh_len);

    return 0;
}
