#ifndef PARTICLE_EFFECTS_H
#define PARTICLE_EFFECTS_H

#include <particle_array.h>

particle_effect newton_step_effect();
particle_effect linear_force_effect(float x, float y, float z);
particle_effect gravitational_force_effect(float x, float y, float z, float mu);
particle_effect plane_bounce_effect(float x, float y, float z, float d);
particle_effect sphere_bounce_effect(float x, float y, float z, float r);

#endif
