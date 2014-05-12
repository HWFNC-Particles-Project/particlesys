#ifndef PARTICLE_ARRAY_H
#define PARTICLE_ARRAY_H

#include <stdint.h>
#include <stddef.h>

/* 
 * the particle type consisting of position, velocity, mass and charge.
 * charge may not have a predefined meaning but mostly serves as an
 * additional simulation parameter and more importantly pads the
 * particle type to 8 floats (half a cache line or two sse vectors or
 * one avx vector).
 */
typedef union particle_t {
    struct {
        float position[3];
        float mass;
        float velocity[3];
        float charge;
    };
    float array[8];
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
 * create an empty particle array. returns 0 on success
 */ 
int particle_array_create(particle_array *array);

/*
 * copies the content of src over to array. array must be created beforehand. returns 0 on success
 */ 
int particle_array_copy(particle_array *array, const particle_array *src);

/*
 * frees a particle array. returns 0 on success
 */ 
int particle_array_destroy(particle_array *array);

/*
 * Compares two particle arrays. Returns 0 when both hold the same particles with tolerance eps.
 */
int particle_array_compare(particle_array *arr1, particle_array *arr2, float eps);

/*
 * returns the current size of the particle array.
 */ 
size_t particle_array_size(const particle_array *array);

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

#endif
