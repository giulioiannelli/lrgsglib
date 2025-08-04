#include "LRGSG_customs.h"
#include "sfmtrng.h"
#include <math.h>
#include <sys/types.h>

#ifndef __LRGSGVMLIB_H_INC__
#define __LRGSGVMLIB_H_INC__

void voter_model_1step(size_t nd, spin_tp s, size_tp nlen, size_tp *neighs, double_p *edgl);
void voter_model_Nstep(size_t N, spin_tp s, size_tp nlen, size_tp *neighs, double_p *edgl);


#endif /* __LRGSGVMLIB_H_INC__ */