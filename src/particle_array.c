#include "particle_array.h"
#include <stdlib.h>

int particle_array_create(particle_array *array) {
    array->size = 0;
    array->capacity = 4;
    array->particles = malloc(array->capacity*sizeof(particle));
    return array->particles == NULL;
}

int particle_array_destroy(particle_array *array) {
    free(array->particles);
    return 0;
}

size_t particle_array_size(particle_array *array) {
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
        if(particle_array_reserve(array, array->capacity/2*3)) {
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
