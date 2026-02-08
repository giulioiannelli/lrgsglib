#include "LRGSG_binsys.h"

#ifndef __LRGSG_CONTDYNSYS_H_INC__
#define __LRGSG_CONTDYNSYS_H_INC__

/* ------------------------------------------------------------------ */
/* State I/O for float64 states                                        */
/* ------------------------------------------------------------------ */
void read_state_f64(const char *fname, size_t N, double *s);
void write_state_f64(FILE *f, size_t N, const double *s);

/* ------------------------------------------------------------------ */
/* ODE integrators (operate on length-N float64 arrays)               */
/* ------------------------------------------------------------------ */

/* Function pointer type for RHS:  rhs(N, s, dsdt, t, params)         */
typedef void (*rhs_func_t)(size_t N, const double *s, double *dsdt,
                           double t, void *params);

void euler_step(size_t N, double *s, rhs_func_t rhs, double t, double dt,
                void *params);
void rk4_step(size_t N, double *s, rhs_func_t rhs, double t, double dt,
              void *params);

/* ------------------------------------------------------------------ */
/* Generic helpers                                                     */
/* ------------------------------------------------------------------ */
void wrap_to_2pi(size_t N, double *s);
double calc_kuramoto_order_param(size_t N, const double *theta);
double calc_mean_concentration(size_t N, const double *u);

#endif /* __LRGSG_CONTDYNSYS_H_INC__ */
