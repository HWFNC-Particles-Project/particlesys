#ifndef PARTICLE_DRAW_H
#define PARTICLE_DRAW_H

#include <particle_array.h>

int particle_vis_init();
void particle_vis_draw(particle_array *particles);
void particle_vis_deinit();

#endif
