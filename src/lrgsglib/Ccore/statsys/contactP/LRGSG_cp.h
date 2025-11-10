#ifndef __LRGSGCPLIB_H_INC__
#define __LRGSGCPLIB_H_INC__

#include "LRGSG_binsys.h"

#define CP_DIR "%s/lrgsg/%s/"
#define SINI_FNAME CP_DIR "s_" PSTR "%s" BINX
#define SOUT_FNAME CP_DIR "sout_" PSTR "%s" BINX

/* Contact-process helper routines */
double cp_infection_rate(size_t node, spin_tp state, size_tp degree, NodeEdges edges);
double cp_recovery_rate(size_t node, double mu, spin_tp state, size_tp degree, NodeEdges edges);

#endif /* __LRGSGCPLIB_H_INC__ */
