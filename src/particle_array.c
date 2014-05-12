#include "particle_array.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

int particle_array_create(particle_array *array) {
    array->size = 0;
    array->capacity = 8;
    array->particles = malloc(array->capacity*sizeof(particle));
    return array->particles == NULL;
}

int particle_array_destroy(particle_array *array) {
    free(array->particles);
    return 0;
}

int particle_array_copy(particle_array *array, const particle_array *src) {
    // check capacity:
    if (src->capacity > array->capacity) {
        // make bigger:
        if (array->particles != NULL) free(array->particles);
        array->capacity = src->capacity;
        array->particles = malloc(array->capacity*sizeof(particle));
        if (array->particles == NULL) return 1;
    }
    // make same size:
    array->size = src->size;
    // copy:
    memcpy(array->particles, src->particles, array->size * sizeof(particle));
    return 0;
}

int particle_array_compare(particle_array *arr1, particle_array *arr2, float eps) {
    size_t i;
    if      (arr1->size > arr2->size) return 1;
    else if (arr1->size < arr2->size) return -1;
    size_t s = arr1->size;
    // compare each particle:
    for (i = 0; i < s; ++i) {
        size_t j;
        // compare each component:
        for (j = 0; j < sizeof(particle) / sizeof(float); ++j) {
            float d = arr1->particles[i].array[j] - arr2->particles[i].array[j];
            // with tolerance eps:
            if (d > eps) {
                fprintf(stderr, "[particle_array_compare] error: particle %d is different (component %d, delta %g)\n", (int)i, (int)j, (double)d);
                return 1;
            } else if (d < -eps) {
                fprintf(stderr, "[particle_array_compare] error: particle %d is different (component %d, delta %g)\n", (int)i, (int)j, (double)d);
                return -1;
            }
        }
    }
    return 0;
}

size_t particle_array_size(const particle_array *array) {
    return array->size;
}

int particle_array_reserve(particle_array *array, size_t capacity) {
    if(capacity < array->capacity) {
        return 0;
    }
    particle *newparticles = realloc(array->particles, capacity*sizeof(particle));
    if(newparticles != NULL) {
        array->capacity = capacity;
        array->particles = newparticles;
        return 0;
    } else {
        return 1;
    }
}

int particle_array_add(particle_array *array, particle p) {
    if(array->size == array->capacity) {
        if(particle_array_reserve(array, array->capacity*2)) {
            return 1;
        }
    }
    array->particles[array->size++] = p;
    return 0;
}

void particle_array_set(particle_array *array, size_t i, particle p) {
    array->particles[i] = p;
}

particle particle_array_get(particle_array *array, size_t i) {
    return array->particles[i];
}
