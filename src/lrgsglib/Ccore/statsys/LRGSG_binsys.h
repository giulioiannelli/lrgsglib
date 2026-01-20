#include "LRGSG_customs.h"
#include "sfmtrng.h"
#include <math.h>
#include <sys/types.h>

#ifndef __LRGSGBINSYSLIB_H_INC__
#define __LRGSGBINSYSLIB_H_INC__

/* Common string format definitions for binary dynamics systems */
#define PSTR "p=%.3g"
#define GRPH_DIR "%s/graph/%s/"
#define EDGL_FNAME GRPH_DIR "edgelist_" PSTR "%s" BINX

#endif /* __LRGSGBINSYSLIB_H_INC__ */
