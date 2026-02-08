#include "LRGSG_binsys.h"

#ifndef __LRGSG_VECDYNSYS_H_INC__
#define __LRGSG_VECDYNSYS_H_INC__

/* ------------------------------------------------------------------ */
/* State I/O for int32 and float64 vector states                       */
/* ------------------------------------------------------------------ */
void read_state_i32(const char *fname, size_t N, int32_t *s);
void write_state_i32(FILE *f, size_t N, const int32_t *s);
void read_state_f64_vec(const char *fname, size_t N, size_t d, double *s);
void write_state_f64_vec(FILE *f, size_t N, size_t d, const double *s);

/* ------------------------------------------------------------------ */
/* Vector operations                                                   */
/* ------------------------------------------------------------------ */
void random_unit_vector_3d(double *v);
void random_rotation_3d(const double *v_in, double *v_out, double delta);
double dot3(const double *a, const double *b);
void normalize3(double *v);

/* ------------------------------------------------------------------ */
/* Discrete state helpers                                              */
/* ------------------------------------------------------------------ */
int kronecker_delta_i32(int32_t a, int32_t b);

#endif /* __LRGSG_VECDYNSYS_H_INC__ */
