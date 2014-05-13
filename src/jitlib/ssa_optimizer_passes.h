#ifndef SSA_OPTIMIZER_PASSES_H
#define SSA_OPTIMIZER_PASSES_H

#include "ssa_def.h"

void ssa_mark_dead(ssa_block *block);
void ssa_remap_duplicates(ssa_block *block, size_t window);
void ssa_fuse_load_store(ssa_block *block);
void ssa_fold_constants(ssa_block *block);
void ssa_schedule(ssa_block *block, const ssa_scheduling_info *info);

#endif
