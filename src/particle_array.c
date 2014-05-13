#include "particle_array.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "malloc_align.h"

int particle_array_realloc(particle_array *array) {
    size_t align = 64;
    // alloc:
    void *mem = NULL;
    particle *particles = (particle *)malloc_align(array->capacity*sizeof(particle), 64, &mem);
    // check if we have to copy anything:
    if (array->mem != NULL) {
        if (mem != NULL) memcpy(particles, array->particles, sizeof(particle) * array->size);
        // free:
        free(array->mem);
        array->mem = NULL;
        array->particles = NULL;
    }
    if (mem != NULL) {
        array->mem       = mem;
        array->particles = particles;
        return 0;
    } else {
        array->mem = NULL;
        array->particles = NULL;
        array->capacity = 0;
        array->size = 0;
        return -1;
    }
}

int particle_array_create(particle_array *array) {
    array->size = 0;
    array->capacity = 8;
    array->mem = NULL;
    array->particles = NULL;
    return particle_array_realloc(array);
}

int particle_array_destroy(particle_array *array) {
    free(array->mem);
    array->mem = NULL;
    array->particles = NULL;
    return 0;
}

int particle_array_copy(particle_array *array, const particle_array *src) {
    // check capacity:
    if (src->capacity > array->capacity) {
        // make bigger:
        array->capacity = src->capacity;
        if (particle_array_realloc(array)) return 1;
    }
    // make same size:
    array->size = src->size;
    // copy:
    memcpy(array->particles, src->particles, array->size * sizeof(particle));
    return 0;
}

int particle_array_compare(particle_array *arr1, particle_array *arr2, float eps, float rel_eps) {
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
            float mean = (arr1->particles[i].array[j] + arr2->particles[i].array[j]) / 2.0;
            float d_rel = d;
            if (mean != 0.0f) {
                d_rel = d / mean;
            }
            // with tolerance eps:
            if (d_rel > rel_eps && d > eps) {
                fprintf(stderr, "[particle_array_compare] error: particle %d is different (component %d, delta %g, rel %g)\n", (int)i, (int)j, (double)d, (double)d_rel);
                return 1;
            } else if (d_rel < -rel_eps && d < -eps) {
                fprintf(stderr, "[particle_array_compare] error: particle %d is different (component %d, delta %g, rel %g)\n", (int)i, (int)j, (double)d, (double)d_rel);
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
    array->capacity = capacity;
    return particle_array_realloc(array);
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
