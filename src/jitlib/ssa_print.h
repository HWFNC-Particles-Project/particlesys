#ifndef SSA_PRINT_H
#define SSA_PRINT_H

#include "ssa_def.h"

void ssa_print_op(ssa_block *block, size_t i);
void ssa_print(ssa_block *block);

#endif
