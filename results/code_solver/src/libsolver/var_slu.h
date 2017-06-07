#pragma once

/**********************************************************************
                  Global variables for SuperLU Poisson solver
 **********************************************************************/ 
#include "use_flags.h"
#include "extern_define_solver.h"

#if !NTRL_ONLY
#include <superlu/slu_ddefs.h>

XTRNSLU SuperMatrix  *L, *U, *B;
XTRNSLU int      *perm_r; /* row permutations from partial pivoting */
XTRNSLU int      *perm_c; /* column permutation vector */
XTRNSLU SuperLUStat_t stat_SLU;
XTRNSLU int info; 
#endif
