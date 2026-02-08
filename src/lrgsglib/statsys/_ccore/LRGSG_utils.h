#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include "LRGSG_customs.h"
#include "sfmtrng.h"

#ifndef __LRSGLIB_H_INC__
#define __LRSGLIB_H_INC__

/* utils.c */

extern void print_cwd(FILE *output_stream);
extern void fill_array_with_value(double *arr, size_t start, size_t end, double value);
extern double* logspace(double start, double stop, int num);
extern int* logspace_int(double stop, int* num);
extern double sum_vs(size_t n, double *v);
extern double sum_vs_2(size_t n, double *v);

extern uint32_t softplus_u32(int32_t x);

extern bool strsme(const char *str1, const char *str2);
extern void strrmac(char *__strdst, const char *__strsrc, const char __chrctr);
extern void strrmuc(char *__strdst, const char *__strsrc, const char __chrctr);
extern bool strbws(char *__strsrc, const char *__strbegin);
extern size_t strtozu(const char *s);
extern uint32_t strtou32(const char *s);

extern bool __feexist(const char *fn);
extern bool __fnexist(const char *fn);
extern void __fopen(FILE **f, const char *fn, const char *md);
extern void __fread_check(size_t __frdval, size_t __frdcnt);
extern void __popen(FILE **p, const char *fn, const char *md);
extern void __fgets_row(FILE **fc, char *row);
extern void prepend(char *s, const char *t);

void __challoc(void *__ptr);
void *__chMalloc(size_t n);
void *__chCalloc(size_t __nmemb, size_t __size);

extern pid_t call(char *argv[]);
extern void __wait_childs(void);

char *rand_string(char *str, size_t size);
extern Edges __read_bin_EdgeList__(const char *filename, size_t *edge_count);
extern void process_edges(const char *filename, size_t N, Edges *edges,
    NodesEdges *node_edges, size_tp *neigh_len);

extern float calculateFractionEqual(int8_t *array1, int8_t *array2, size_t size);

extern void __make_adj_from_tmp(size_t i, size_t j, double tmp, double_p **adj);
extern void __fill_adj__(FILE **f, size_t N, double_p **adj);
extern void __fill_edgl_read__(FILE **f, size_t N, double_p **edgl, 
    size_tp **neighs, size_tp *neigh_len);
extern void __fill_edgl_make__(FILE **f, size_t N, double_p **adj, 
    double_p **edgl, size_tp **neighs, size_tp *neigh_len);

extern void build_str_id(const char *_str_id,
                 char *str_id,
                 size_t str_sz);


double neigh_weight_magn(size_t nd, size_t n_nn, spin_tp s, size_tp *neighs,
                         double_p *edgl);
double neighWeight_magn(NodeEdges ne, spin_tp s, size_t n_nn);
double calc_magn(size_t N, spin_tp s);
double calc_ext_magn(size_t N, spin_tp s);
double calc_ext_magn2(size_t N, spin_tp s);
void flip_spin(size_t nd, spin_tp s);
double calc_clust_magn(size_t cli_l, size_tp cli, spin_tp s);
double calc_totEnergy(size_t N, spin_tp s, size_tp nlen, NodesEdges ne);
double calc_ising_energy_normalized(size_t N, spin_tp s, size_tp nlen, NodesEdges ne);
double calc_ising_energy_fast(size_t N, spin_tp s, size_tp nlen, NodesEdges ne);
bool glauber_isStableAtZeroTemp(size_t N, spin_tp s, size_tp nlen, 
    NodesEdges ne);

#endif /* __LRSGLIB_H_INC__ */