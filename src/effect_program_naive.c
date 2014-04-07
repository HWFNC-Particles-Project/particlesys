#include "effect_program_naive.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*
 * particle effects consists of a function (apply) and userdata for
 * effect specific data that is passed as second function to apply
 */
typedef struct particle_effect_naive_t {
    void (*apply)(particle *, void *, float);
    void *userdata;
} particle_effect_naive;


void effect_program_naive_compile(      effect_program *self, const effect_desc *desc);
void effect_program_naive_execute(const effect_program *self,    particle_array *arr, float dt);
void effect_program_naive_destroy(effect_program *self);

particle_effect_naive linear_accel_effect(float x, float y, float z);
particle_effect_naive linear_force_effect(float x, float y, float z);
particle_effect_naive gravitational_force_effect(float x, float y, float z, float mu);
particle_effect_naive plane_bounce_effect(float x, float y, float z, float d, float a);
particle_effect_naive sphere_bounce_effect(float x, float y, float z, float r, float a);
particle_effect_naive newton_step_effect();

void linear_accel_apply(particle *p, void *data0, float dt);
void linear_force_apply(particle *p, void *data0, float dt);
void gravitational_force_apply(particle *p, void *data0, float dt);
void plane_bounce_apply(particle *p, void *data0, float dt);
void sphere_bounce_apply(particle *p, void *data0, float dt);
void newton_step_apply(particle *p, void *data0, float dt);



int effect_program_create_naive(effect_program *p) {
    memset(p, 0, sizeof(effect_program));
    p->compile = effect_program_naive_compile;
    p->execute = effect_program_naive_execute;
    p->destroy = effect_program_naive_destroy;
    p->usr = NULL;
    return 0;
}

void effect_program_naive_destroy(effect_program *self) {
    if (self->usr != NULL) {
        // free last compilation result:
        particle_effect_naive *effects = (particle_effect_naive *)self->usr;
        for(size_t j = 0; effects[j].apply != NULL; ++j) {
            if (effects[j].userdata != NULL) {
                free(effects[j].userdata);
            }
        }
        free(self->usr);
        self->usr = NULL;
    }
}

void effect_program_naive_compile(effect_program *self, const effect_desc *desc) {
    if (self->usr != NULL) {
        effect_program_naive_destroy(self);
    }
    // compile ...
    size_t count = desc->size;
    self->usr = malloc((count + 1) * sizeof(particle_effect_naive));
    particle_effect_naive *effects = (particle_effect_naive *)self->usr;
    for(size_t i = 0; i < count; ++i) {
        const effect_desc_ele *el = &desc->elements[i];
        switch(el->type) {
            case EFFECT_TYPE_LINEAR_ACCEL:
                effects[i] = linear_accel_effect(       el->float_usr[0],
                                                        el->float_usr[1],
                                                        el->float_usr[2]);
                break;
            case EFFECT_TYPE_LINEAR_FORCE:
                effects[i] = linear_force_effect(       el->float_usr[0],
                                                        el->float_usr[1],
                                                        el->float_usr[2]);
                break;
            case EFFECT_TYPE_CENTRAL_FORCE:
                effects[i] = gravitational_force_effect(el->float_usr[0],
                                                        el->float_usr[1],
                                                        el->float_usr[2],
                                                        el->float_usr[3]);
                break;
            case EFFECT_TYPE_PLANE_BOUNCE:
                effects[i] = plane_bounce_effect(       el->float_usr[0],
                                                        el->float_usr[1],
                                                        el->float_usr[2],
                                                        el->float_usr[3],
                                                        el->float_usr[4]);
                break;
            case EFFECT_TYPE_SPHERE_BOUNCE:
                effects[i] = sphere_bounce_effect(      el->float_usr[0],
                                                        el->float_usr[1],
                                                        el->float_usr[2],
                                                        el->float_usr[3],
                                                        el->float_usr[4]);
                break;
            case EFFECT_TYPE_NEWTON_STEP:
                effects[i] = newton_step_effect();
                break;
        }
    }
    effects[count].apply = NULL;
    effects[count].userdata = NULL;

}

void effect_program_naive_execute(const effect_program *self, particle_array *arr, float dt) {
    if (self->usr == NULL) {
        // fail, no compilation result.
        return;
    }
    particle_effect_naive *effects = (particle_effect_naive *)self->usr;
    for(size_t i = 0; i < arr->size; ++i) {
        for(size_t j = 0; effects[j].apply != NULL; ++j) {
            effects[j].apply(&arr->particles[i], effects[j].userdata, dt);
        }
    }
}

void linear_accel_apply(particle *p, void *data0, float dt) {
    float *data = data0;
    p->velocity[0] += dt*data[0];
    p->velocity[1] += dt*data[1];
    p->velocity[2] += dt*data[2];
}

particle_effect_naive linear_accel_effect(float x, float y, float z) {
    particle_effect_naive result;
    result.apply = linear_accel_apply;
    float *data = malloc(3*sizeof(float));
    data[0] = x;
    data[1] = y;
    data[2] = z;
    result.userdata = data;
    return result;
}

void linear_force_apply(particle *p, void *data0, float dt) {
    float *data = data0;
    p->velocity[0] += dt*data[0] / p->mass;
    p->velocity[1] += dt*data[1] / p->mass;
    p->velocity[2] += dt*data[2] / p->mass;
}

particle_effect_naive linear_force_effect(float x, float y, float z) {
    particle_effect_naive result;
    result.apply = linear_force_apply;
    float *data = malloc(3*sizeof(float));
    data[0] = x;
    data[1] = y;
    data[2] = z;
    result.userdata = data;
    return result;
}

void gravitational_force_apply(particle *p, void *data0, float dt) {
    float *data = data0;
    float diff[3];
    float mu = data[3];
    diff[0] = p->position[0]-data[0];
    diff[1] = p->position[1]-data[1];
    diff[2] = p->position[2]-data[2];
    float r = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
    float r3 = 1.0f/(r*r*r);
    p->velocity[0] += dt*mu*diff[0]*r3;
    p->velocity[1] += dt*mu*diff[1]*r3;
    p->velocity[2] += dt*mu*diff[2]*r3;
}

particle_effect_naive gravitational_force_effect(float x, float y, float z, float mu) {
    particle_effect_naive result;
    result.apply = gravitational_force_apply;
    float *data = malloc(4*sizeof(float));
    data[0] = x;
    data[1] = y;
    data[2] = z;
    data[3] = mu;
    result.userdata = data;
    return result;
}

void plane_bounce_apply(particle *p, void *data0, float dt) {
    (void) dt;
    float *data = data0;
    float dist = data[0]*p->position[0] + data[1]*p->position[1] + data[2]*p->position[2];
    float vnormal = data[0]*p->velocity[0] + data[1]*p->velocity[1] + data[2]*p->velocity[2];
    float d = data[3];
    if(dist<d && vnormal<0) {
        p->velocity[0] -= (1.0+data[4])*vnormal*data[0];
        p->velocity[1] -= (1.0+data[4])*vnormal*data[1];
        p->velocity[2] -= (1.0+data[4])*vnormal*data[2];
    }
}

particle_effect_naive plane_bounce_effect(float x, float y, float z, float d, float a) {
    particle_effect_naive result;
    result.apply = plane_bounce_apply;
    float *data = malloc(5*sizeof(float));
    float r = sqrtf(x*x + y*y + z*z);
    data[0] = x/r;
    data[1] = y/r;
    data[2] = z/r;
    data[3] = d/r;
    data[4] = a;
    result.userdata = data;
    return result;
}

void sphere_bounce_apply(particle *p, void *data0, float dt) {
    (void) dt;
    float *data = data0;
    float normal[3];
    normal[0] = p->position[0]-data[0];
    normal[1] = p->position[1]-data[1];
    normal[2] = p->position[2]-data[2];
    double r = sqrtf(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
    normal[0] /= r;
    normal[1] /= r;
    normal[2] /= r;
    float vnormal = normal[0]*p->velocity[0] + normal[1]*p->velocity[1] + normal[2]*p->velocity[2];
    float d = data[3];
    if(r<d && vnormal<0) {
        p->velocity[0] -= (1.0+data[4])*vnormal*normal[0];
        p->velocity[1] -= (1.0+data[4])*vnormal*normal[1];
        p->velocity[2] -= (1.0+data[4])*vnormal*normal[2];
    }
}

particle_effect_naive sphere_bounce_effect(float x, float y, float z, float r, float a) {
    particle_effect_naive result;
    result.apply = sphere_bounce_apply;
    float *data = malloc(5*sizeof(float));
    data[0] = x;
    data[1] = y;
    data[2] = z;
    data[3] = r;
    data[4] = a;
    result.userdata = data;
    return result;
}

void newton_step_apply(particle *p, void *data0, float dt) {
    (void) data0;
    p->position[0] += dt*p->velocity[0];
    p->position[1] += dt*p->velocity[1];
    p->position[2] += dt*p->velocity[2];
}

particle_effect_naive newton_step_effect() {
    particle_effect_naive result;
    result.apply = newton_step_apply;
    result.userdata = NULL;
    return result;
}

