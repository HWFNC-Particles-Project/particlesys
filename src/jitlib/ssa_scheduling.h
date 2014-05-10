#ifndef SSA_SCHEDULING_H
#define SSA_SCHEDULING_H

#include "ssa_def.h"

static ssa_scheduling_info sandybridge[] = {
    [SSA_NONE]   = {0 ,0 ,0},
    [SSA_ADD]    = {3 ,1 ,1},
    [SSA_MUL]    = {5 ,1 ,0},
    [SSA_MIN]    = {3 ,1 ,0},
    [SSA_MAX]    = {3 ,1 ,0},
    [SSA_CMPEQ]  = {3, 1, 0},
    [SSA_CMPLE]  = {3, 1, 0},
    [SSA_CMPLT]  = {3, 1, 0},
    [SSA_CMPNEQ] = {3, 1, 0},
    [SSA_CMPNLE] = {3, 1, 0},
    [SSA_CMPNLT] = {3, 1, 0},
    [SSA_AND]    = {1, 1, 0},
    [SSA_ANDN]   = {1, 1, 0},
    [SSA_XOR]    = {1, 1, 0},
    [SSA_OR]     = {1, 1, 0},
    [SSA_STORE]  = {0 ,0 ,0},
    [SSA_DIV]    = {14,14,0},
    [SSA_SUB]    = {3 ,1 ,1},
    [SSA_LOAD]   = {0 ,0 ,0},
    [SSA_CONST]  = {1 ,1 ,0},
    [SSA_SQRT]   = {14,14,0},
    [SSA_RSQRT]  = {5 ,1 ,0},
    [SSA_RCP]    = {5 ,1 ,0},
    [SSA_BLEND]  = {2 ,1 ,0},
};

static ssa_scheduling_info ivybridge[] = {
    [SSA_NONE]   = {0 ,0 ,0},
    [SSA_ADD]    = {3 ,1 ,1},
    [SSA_MUL]    = {5 ,1 ,0},
    [SSA_MIN]    = {3 ,1 ,0},
    [SSA_MAX]    = {3 ,1 ,0},
    [SSA_CMPEQ]  = {3, 1, 0},
    [SSA_CMPLE]  = {3, 1, 0},
    [SSA_CMPLT]  = {3, 1, 0},
    [SSA_CMPNEQ] = {3, 1, 0},
    [SSA_CMPNLE] = {3, 1, 0},
    [SSA_CMPNLT] = {3, 1, 0},
    [SSA_AND]    = {1, 1, 0},
    [SSA_ANDN]   = {1, 1, 0},
    [SSA_XOR]    = {1, 1, 0},
    [SSA_OR]     = {1, 1, 0},
    [SSA_STORE]  = {0 ,0 ,0},
    [SSA_DIV]    = {13,7 ,0},
    [SSA_SUB]    = {3 ,1 ,1},
    [SSA_LOAD]   = {0 ,0 ,0},
    [SSA_CONST]  = {1 ,1 ,0},
    [SSA_SQRT]   = {13,7 ,0},
    [SSA_RSQRT]  = {5 ,1 ,0},
    [SSA_RCP]    = {5 ,1 ,0},
    [SSA_BLEND]  = {2 ,1 ,0},
};

#endif
