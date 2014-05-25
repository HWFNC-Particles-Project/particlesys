#ifndef __effect_program_jit_H
#define __effect_program_jit_H

#include "effect_program.h"
#include "particle_array.h"

#define JIT_AVX1 0x00
#define JIT_AVX4 0x01
#define JIT_AVX8 0x02


#define JIT_O0 0x00
#define JIT_O1 0x10
#define JIT_O2 0x20
#define JIT_O3 0x30

int effect_program_create_jit(effect_program *p, int level);

#endif  // __effect_program_jit_H
