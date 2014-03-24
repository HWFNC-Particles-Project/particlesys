#ifndef PARTICLE_ARRAY_H
#define PARTICLE_ARRAY_H

#include <stdlib.h>

/* 
 * the particle type consisting of position, velocity, mass and charge.
 * charge may not have a predefined meaning but mostly serves as an
 * additional simulation parameter and more importantly pads the
 * particle type to 8 floats (half a cache line or two sse vectors or
 * one avx vector).
 */
typedef struct particle_t {
    float position[3];
    float mass;
    float velocity[3];
    float charge;
} particle;

/*
 * dynamically sized particle array that supports bulk application of
 * effects
 */
typedef struct particle_array_t {
    size_t size, capacity;
    particle *particles;
} particle_array;

/*
 * particle effects consiste of a function (apply) and userdata for
 * effect specific data that is passed as second function to apply
 */ 
typedef struct particle_effect_t {
    void (*apply)(particle*, void*, float);
    void *userdata;
} particle_effect;

/*
 * frees a zero terminated effect array by calling free on all
 * the userdata pointers
 */ 
void particle_effect_free(particle_effect *effects);

/*
 * create an empty particle array. returns 0 on success
 */ 
int particle_array_create(particle_array *array);

/*
 * frees a particle array. returns 0 on success
 */ 
int particle_array_destroy(particle_array *array);

/*
 * returns the current size of the particle array.
 */ 
size_t particle_array_size(particle_array *array);

/*
 * ensures a minimum capacity. returns 0 on success.
 */ 
int particle_array_reserve(particle_array *array, size_t capacity);

/*
 * add a particle to the particle array and resizes if required.
 * returns 0 on success.
 */ 
int particle_array_add(particle_array *array, particle p);

/*
 * set the value of particle at positions i to p.
 */ 
void particle_array_set(particle_array *array, size_t i, particle p);

/*
 * returns the particle at positions i to p.
 */ 
particle particle_array_get(particle_array *array, size_t i);

/*
 * applies a zero terminated array of effects to the array.
 */ 
void particle_array_apply_effects(particle_array *array, particle_effect *effects, float dt);

#endif
