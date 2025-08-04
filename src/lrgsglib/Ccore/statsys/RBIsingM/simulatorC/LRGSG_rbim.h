#include "LRGSG_customs.h"
#include "sfmtrng.h"
#include <math.h>
#include <sys/types.h>
//
#ifndef __LRGSGRBIMLIB_H_INC__
#define __LRGSGRBIMLIB_H_INC__
//
#define BOLTZMANN_FACTOR(DE, T) exp(-DE / T)
#define T_THERM_STEP  (size_t)(thSTEP * N)
#define T_EQ_STEP (size_t)(eqSTEP * N)
#define T_STEPS (T_THERM_STEP + T_EQ_STEP)
//
#define PSTR "p=%.3g"
#define TSTR "T=%.3g"
#define PSTR_TSTR PSTR "_" TSTR
//
#define ISNG_DIR "%s/ising/%s/"
#define SINI_FNAME ISNG_DIR     "s_"        PSTR    "%s" BINX
#define CLID_FNAME ISNG_DIR     "cl%zu_"    PSTR    "%s" BINX
#define HFIELD_FNAME ISNG_DIR   "h_"        PSTR    "%s" BINX
//
#define CLOUT_FNAME ISNG_DIR    "outcl%zu_" PSTR_TSTR   "%s" "%s"
#define ENE_FNAME ISNG_DIR      "ene_"      PSTR_TSTR   "%s" BINX
#define ENEF_FNAME ISNG_DIR     "ene_"      PSTR_TSTR   "%s" TXTX
#define SOUT_FNAME ISNG_DIR     "sout_"     PSTR_TSTR   "%s" BINX
#define MAGN_FNAME ISNG_DIR     "m_"        PSTR_TSTR   "%s" BINX
#define MAGNF_FNAME ISNG_DIR    "m_"        PSTR_TSTR   "%s" TXTX
//
#define GRPH_DIR "%s/graph/%s/"
#define EDGL_FNAME GRPH_DIR "edgelist_" PSTR "%s" BINX
#define EIGV_FNAME GRPH_DIR "eigV%d_"   PSTR "%s" BINX
#define ADJ_FNAME GRPH_DIR "adj_"       PSTR "%s" BINX
/* * Glauber–Metropolis kernels (T > 0)
 *  ----------------------------------------------------
 *
 *  g_h == NULL      →   no field
 *  g_h != NULL      →   site–dependent field h[i]
 * */
typedef void (*glauMetroFunc)(size_t, double, spin_tp, size_t, NodeEdges);
//
extern sfmt_t sfmt;
extern uint32_t *seed_rand;
extern double thSTEP;
extern double eqSTEP;
//
void initialize_glauberMetropolis(double T, const double *h);
void glauberMetropolis_1step(size_t nd, double T, spin_tp s, size_t nlen,
    NodeEdges ne);
void glauberMetropolis_1step_T0(size_t nd, double T, spin_tp s, size_t nlen,
    NodeEdges ne);
void glauberMetropolis_Nstep(size_t N, double T, spin_tp s, size_tp nlen,
    NodesEdges ne, const char *update_mode);
//
void glauber_metropolis_1step(size_t nd, double T, spin_tp s, size_tp nlen,
                              size_tp *neighs, double_p *edgl);
void glauber_metropolis_Nstep(size_t N, double T, spin_tp s, size_tp nlen,
                              size_tp *neighs, double_p *edgl);
void glauber_metropolis_Nstep_mode(size_t N, double T, spin_tp s, size_tp nlen,
                       size_tp *neighs, double_p *edgl, const char *update_mode);
double calc_energy_full(size_t N, spin_tp s, size_tp nlen, size_tp *neighs,
                        double_p *edgl);

#endif /* __LRGSGRBIMLIB_H_INC__ */