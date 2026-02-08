#include "LRGSG_vecdynsys.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"
#include <math.h>

/* ================================================================== */
/* State I/O                                                           */
/* ================================================================== */

void read_state_i32(const char *fname, size_t N, int32_t *s) {
    FILE *f;
    __fopen(&f, fname, "rb");
    __fread_check(fread(s, sizeof(int32_t), N, f), N);
    fclose(f);
}

void write_state_i32(FILE *f, size_t N, const int32_t *s) {
    fwrite(s, sizeof(int32_t), N, f);
    fflush(f);
}

void read_state_f64_vec(const char *fname, size_t N, size_t d, double *s) {
    FILE *f;
    __fopen(&f, fname, "rb");
    __fread_check(fread(s, sizeof(double), N * d, f), N * d);
    fclose(f);
}

void write_state_f64_vec(FILE *f, size_t N, size_t d, const double *s) {
    fwrite(s, sizeof(double), N * d, f);
    fflush(f);
}

/* ================================================================== */
/* Vector operations                                                   */
/* ================================================================== */

void random_unit_vector_3d(double *v) {
    /* Marsaglia method for uniform points on S^2 */
    double x1, x2, s_sq;
    do {
        x1 = 2.0 * genrand_real2() - 1.0;
        x2 = 2.0 * genrand_real2() - 1.0;
        s_sq = x1 * x1 + x2 * x2;
    } while (s_sq >= 1.0);
    double factor = 2.0 * sqrt(1.0 - s_sq);
    v[0] = x1 * factor;
    v[1] = x2 * factor;
    v[2] = 1.0 - 2.0 * s_sq;
}

void random_rotation_3d(const double *v_in, double *v_out, double delta) {
    /* Small random rotation of v_in by angle up to delta */
    double perturb[3];
    for (int k = 0; k < 3; ++k)
        perturb[k] = v_in[k] + delta * (2.0 * genrand_real2() - 1.0);
    normalize3(perturb);
    v_out[0] = perturb[0];
    v_out[1] = perturb[1];
    v_out[2] = perturb[2];
}

double dot3(const double *a, const double *b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void normalize3(double *v) {
    double norm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (norm > 0.0) {
        v[0] /= norm;
        v[1] /= norm;
        v[2] /= norm;
    }
}

/* ================================================================== */
/* Discrete state helpers                                              */
/* ================================================================== */

int kronecker_delta_i32(int32_t a, int32_t b) {
    return (a == b) ? 1 : 0;
}
