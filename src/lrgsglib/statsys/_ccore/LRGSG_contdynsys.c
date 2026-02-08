#include "LRGSG_contdynsys.h"
#include "LRGSG_utils.h"
#include <string.h>
#include <math.h>

/* ================================================================== */
/* State I/O                                                           */
/* ================================================================== */

void read_state_f64(const char *fname, size_t N, double *s) {
    FILE *f;
    __fopen(&f, fname, "rb");
    __fread_check(fread(s, sizeof(double), N, f), N);
    fclose(f);
}

void write_state_f64(FILE *f, size_t N, const double *s) {
    fwrite(s, sizeof(double), N, f);
    fflush(f);
}

/* ================================================================== */
/* ODE integrators                                                     */
/* ================================================================== */

void euler_step(size_t N, double *s, rhs_func_t rhs, double t, double dt,
                void *params) {
    double *dsdt = __chMalloc(N * sizeof(double));
    rhs(N, s, dsdt, t, params);
    for (size_t i = 0; i < N; ++i)
        s[i] += dt * dsdt[i];
    free(dsdt);
}

void rk4_step(size_t N, double *s, rhs_func_t rhs, double t, double dt,
              void *params) {
    double *k1   = __chMalloc(N * sizeof(double));
    double *k2   = __chMalloc(N * sizeof(double));
    double *k3   = __chMalloc(N * sizeof(double));
    double *k4   = __chMalloc(N * sizeof(double));
    double *tmp  = __chMalloc(N * sizeof(double));

    /* k1 */
    rhs(N, s, k1, t, params);

    /* k2 */
    for (size_t i = 0; i < N; ++i) tmp[i] = s[i] + 0.5 * dt * k1[i];
    rhs(N, tmp, k2, t + 0.5 * dt, params);

    /* k3 */
    for (size_t i = 0; i < N; ++i) tmp[i] = s[i] + 0.5 * dt * k2[i];
    rhs(N, tmp, k3, t + 0.5 * dt, params);

    /* k4 */
    for (size_t i = 0; i < N; ++i) tmp[i] = s[i] + dt * k3[i];
    rhs(N, tmp, k4, t + dt, params);

    /* combine */
    for (size_t i = 0; i < N; ++i)
        s[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);

    free(k1); free(k2); free(k3); free(k4); free(tmp);
}

/* ================================================================== */
/* Generic helpers                                                     */
/* ================================================================== */

void wrap_to_2pi(size_t N, double *s) {
    for (size_t i = 0; i < N; ++i) {
        s[i] = fmod(s[i], 2.0 * M_PI);
        if (s[i] < 0.0) s[i] += 2.0 * M_PI;
    }
}

double calc_kuramoto_order_param(size_t N, const double *theta) {
    double sum_cos = 0.0, sum_sin = 0.0;
    for (size_t i = 0; i < N; ++i) {
        sum_cos += cos(theta[i]);
        sum_sin += sin(theta[i]);
    }
    sum_cos /= (double)N;
    sum_sin /= (double)N;
    return sqrt(sum_cos * sum_cos + sum_sin * sum_sin);
}

double calc_mean_concentration(size_t N, const double *u) {
    double sum = 0.0;
    for (size_t i = 0; i < N; ++i)
        sum += u[i];
    return sum / (double)N;
}
