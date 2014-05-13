#include "effect_program_naive.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*
 * particle effects consists of a function (apply) and userdata for
 * effect specific data that is passed as second function to apply
 */
typedef struct particle_effect_naive_t {
    int particles;
    union {
        void (*one)(particle *, void *, float);
        void (*two)(particle *, particle *, void *, float);
    } apply;
    union {
        void (*one)(const particle *, void *, float, performance_count *out);
        void (*two)(const particle *, const particle *, void *, float, performance_count *out);
    } perf_c;
    void *userdata;
} particle_effect_naive;


void effect_program_naive_compile(effect_program *self, const effect_desc *desc);
void effect_program_naive_execute(const effect_program *self, particle_array *arr, float dt);
void effect_program_naive_perf_c(const effect_program *self, const particle_array *arr, float dt, performance_count *out);
void effect_program_naive_destroy(effect_program *self);

static particle_effect_naive linear_accel_effect(float x, float y, float z);
static particle_effect_naive linear_force_effect(float x, float y, float z);
static particle_effect_naive gravitational_force_effect(float x, float y, float z, float mu);
static particle_effect_naive plane_bounce_effect(float x, float y, float z, float d, float a);
static particle_effect_naive sphere_bounce_effect(float x, float y, float z, float r, float a);
static particle_effect_naive pairwise_gravitational_force_effect(float mu);
static particle_effect_naive pairwise_sphere_collision_effect(float radius, float restitution);
static particle_effect_naive newton_step_effect();

int effect_program_create_naive(effect_program *p) {
    memset(p, 0, sizeof(effect_program));
    p->compile = effect_program_naive_compile;
    p->execute = effect_program_naive_execute;
    p->perf_c  = effect_program_naive_perf_c;
    p->destroy = effect_program_naive_destroy;
    p->usr = NULL;
    return 0;
}

void effect_program_naive_destroy(effect_program *self) {
    if (self->usr != NULL) {
        // free last compilation result:
        particle_effect_naive *effects = (particle_effect_naive *)self->usr;
        for(size_t j = 0; effects[j].apply.one != NULL; ++j) {
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

            case EFFECT_TYPE_GRAVITY_FORCE:
                effects[i] = pairwise_gravitational_force_effect(el->float_usr[0]);
                break;

            case EFFECT_TYPE_SPHERE_COLLISION:
                effects[i] = pairwise_sphere_collision_effect(  el->float_usr[0],
                                                                el->float_usr[1]);
                break;
            default:
                break; // unhandled effect
        }
    }
    effects[count].apply.one = NULL;
    effects[count].userdata = NULL;

}

void effect_program_naive_execute(const effect_program *self, particle_array *arr, float dt) {
    if (self->usr == NULL) {
        // fail, no compilation result.
        return;
    }
    particle_effect_naive *effects = (particle_effect_naive *)self->usr;
    for(size_t i = 0; i < arr->size; ++i) {
        for(size_t j = 0; effects[j].apply.one != NULL; ++j) {
            if(effects[j].particles == 1) {
                effects[j].apply.one(&arr->particles[i], effects[j].userdata, dt);
            } else if(effects[j].particles == 2) {
                for(size_t k = 0; k<i; ++k) {
                    effects[j].apply.two(&arr->particles[i], &arr->particles[k], effects[j].userdata, dt);
                }
            }
        }
    }
}

void effect_program_naive_perf_c(const effect_program *self, const particle_array *arr, float dt, performance_count *out) {
    // set all counters to zero:
    memset(out, 0, sizeof(*out));
    if (self->usr == NULL) {
        // fail, no compilation result.
        return;
    }
    particle_effect_naive *effects = (particle_effect_naive *)self->usr;
    for(size_t i = 0; i < arr->size; ++i) {
        for(size_t j = 0; effects[j].apply.one != NULL; ++j) {
            if(effects[j].particles == 1) {
                if (effects[j].perf_c.one != NULL) 
                    effects[j].perf_c.one(&arr->particles[i], effects[j].userdata, dt, out);
            } else if(effects[j].particles == 2) {
                for(size_t k = 0; k<i; ++k) {
                    if (effects[j].perf_c.two != NULL) 
                        effects[j].perf_c.two(&arr->particles[i], &arr->particles[k], effects[j].userdata, dt, out);
                }
            }
        }
    }
}

static void linear_accel_apply(particle *p, void *data0, float dt) {
    float *data = data0;
    p->velocity[0] += dt*data[0];
    p->velocity[1] += dt*data[1];
    p->velocity[2] += dt*data[2];
}

static void linear_accel_perf_c(const particle *p, void *data0, float dt, performance_count *out) {
    (void) p; (void) data0; (void) dt;
    out->add += 3;
    out->mul += 3;
    out->loads += 6;
    out->stores += 3;
}

static particle_effect_naive linear_accel_effect(float x, float y, float z) {
    particle_effect_naive result;
    result.particles = 1;
    result.apply.one  = linear_accel_apply;
    result.perf_c.one = linear_accel_perf_c;
    float *data = malloc(3*sizeof(float));
    data[0] = x;
    data[1] = y;
    data[2] = z;
    result.userdata = data;
    return result;
}

static void linear_force_apply(particle *p, void *data0, float dt) {
    float *data = data0;
    p->velocity[0] += data[0] * (dt/p->mass);
    p->velocity[1] += data[1] * (dt/p->mass);
    p->velocity[2] += data[2] * (dt/p->mass);
}

static void linear_force_perf_c(const particle *p, void *data0, float dt, performance_count *out) {
    (void) p; (void) data0; (void) dt;
    out->add += 3;
    out->mul += 3;
    out->div += 1;  // is probably optimized by compiler.
    out->loads += 7;
    out->stores += 3;
}

static particle_effect_naive linear_force_effect(float x, float y, float z) {
    particle_effect_naive result;
    result.particles = 1;
    result.apply.one  = linear_force_apply;
    result.perf_c.one = linear_force_perf_c;
    float *data = malloc(3*sizeof(float));
    data[0] = x;
    data[1] = y;
    data[2] = z;
    result.userdata = data;
    return result;
}

static void gravitational_force_apply(particle *p, void *data0, float dt) {
    float *data = data0;
    float diff[3];
    float mu = data[3];
    diff[0] = p->position[0]-data[0];
    diff[1] = p->position[1]-data[1];
    diff[2] = p->position[2]-data[2];
    float r = sqrtf(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
    float r3 = 1.0f/(r*r*r);
    p->velocity[0] += dt*mu*diff[0]*r3;
    p->velocity[1] += dt*mu*diff[1]*r3;
    p->velocity[2] += dt*mu*diff[2]*r3;
}

static void gravitational_force_perf_c(const particle *p, void *data0, float dt, performance_count *out) {
    (void) p; (void) data0; (void) dt;
    out->add += 8;
    out->mul += 12; // compiler might optimize some multiplications away (dt * mu)
    out->div += 1;
    out->sqrt += 1;
    out->loads += 10;
    out->stores += 3;
}

static particle_effect_naive gravitational_force_effect(float x, float y, float z, float mu) {
    particle_effect_naive result;
    result.particles = 1;
    result.apply.one  = gravitational_force_apply;
    result.perf_c.one = gravitational_force_perf_c;
    float *data = malloc(4*sizeof(float));
    data[0] = x;
    data[1] = y;
    data[2] = z;
    data[3] = mu;
    result.userdata = data;
    return result;
}

static void plane_bounce_apply(particle *p, void *data0, float dt) {
    (void) dt;
    float *data = data0;
    float dist =    data[0]*p->position[0] + data[1]*p->position[1] + data[2]*p->position[2];
    float vnormal = data[0]*p->velocity[0] + data[1]*p->velocity[1] + data[2]*p->velocity[2];
    float d = data[3];
    if(0.0f<d-dist && vnormal<0.0f) {
        // we are behind plane and velocity is away from the back of the plane
        // invert normal components of velocity:
        p->velocity[0] -= 2.0f*vnormal*data[0];
        p->velocity[1] -= 2.0f*vnormal*data[1];
        p->velocity[2] -= 2.0f*vnormal*data[2];
        // apply decay:
        p->velocity[0] *= data[4];
        p->velocity[1] *= data[4];
        p->velocity[2] *= data[4];
        // move particle to the surface of the plane (works only when vector data[0..2] is normalized.
        p->position[0] += (d-dist)*data[0];
        p->position[1] += (d-dist)*data[1];
        p->position[2] += (d-dist)*data[2];
    }
}

static void plane_bounce_perf_c(const particle *p, void *data0, float dt, performance_count *out) {
    (void) dt;
    float *data = data0;
    float dist =    data[0]*p->position[0] + data[1]*p->position[1] + data[2]*p->position[2];
    float vnormal = data[0]*p->velocity[0] + data[1]*p->velocity[1] + data[2]*p->velocity[2];
    float d = data[3];
    out->add += 4;
    out->mul += 6;
    out->loads += 10;
    out->cmp += 1;
    if(dist<d) {
        out->cmp += 1;
        if (vnormal<0.0f) {
            // branch is taken:
            out->add += 7;
            out->mul += 10;
            out->loads += 1;
            out->stores += 6;
        }
    } else {
        // branch not taken.
    }
}

static particle_effect_naive plane_bounce_effect(float x, float y, float z, float d, float a) {
    particle_effect_naive result;
    result.particles = 1;
    result.apply.one =  plane_bounce_apply;
    result.perf_c.one = plane_bounce_perf_c;
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

static void sphere_bounce_apply(particle *p, void *data0, float dt) {
    (void) dt;
    float *data = data0;
    float normal[3];
    normal[0] = p->position[0]-data[0];
    normal[1] = p->position[1]-data[1];
    normal[2] = p->position[2]-data[2];
    float r = sqrtf(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
    normal[0] /= r;
    normal[1] /= r;
    normal[2] /= r;
    float vnormal = normal[0]*p->velocity[0] + normal[1]*p->velocity[1] + normal[2]*p->velocity[2];
    float d = data[3];
    if(r<d && vnormal<0.0f) {
        // inside sphere and going further inside sphere:
        p->velocity[0] -= 2.0f*vnormal*normal[0];
        p->velocity[1] -= 2.0f*vnormal*normal[1];
        p->velocity[2] -= 2.0f*vnormal*normal[2];

        p->velocity[0] *= data[4];
        p->velocity[1] *= data[4];
        p->velocity[2] *= data[4];

        p->position[0] += (d-r)*normal[0];
        p->position[1] += (d-r)*normal[1];
        p->position[2] += (d-r)*normal[2];
    }
}

static void sphere_bounce_perf_c(const particle *p, void *data0, float dt, performance_count *out) {
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
    out->add += 7;
    out->mul += 6;
    out->div += 3;
    out->sqrt += 1;
    out->loads += 10;
    out->cmp += 1;
    if(r<d) {
        out->cmp += 1;
        if (vnormal<0.0f) {
            // branch is taken:
            out->add += 7;
            out->mul += 10;
            out->loads += 1;
            out->stores += 6;
        }
    } else {
        // branch not taken.
    }
}

static particle_effect_naive sphere_bounce_effect(float x, float y, float z, float r, float a) {
    particle_effect_naive result;
    result.particles = 1;
    result.apply.one =  sphere_bounce_apply;
    result.perf_c.one = sphere_bounce_perf_c;
    float *data = malloc(5*sizeof(float));
    data[0] = x;
    data[1] = y;
    data[2] = z;
    data[3] = r;
    data[4] = a;
    result.userdata = data;
    return result;
}

static void pairwise_gravitational_force_apply(particle *p1, particle *p2, void *data0, float dt) {
    float *data = data0;
    float mu = data[0];
    float diff[3];
    float m1 = p1->mass;
    float m2 = p2->mass;
    mu *= m1*m2;
    diff[0] = p2->position[0]-p1->position[0];
    diff[1] = p2->position[1]-p1->position[1];
    diff[2] = p2->position[2]-p1->position[2];
    float r = sqrtf(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
    // is r^3 because we need to normalize diff.
    mu = mu/(r*r*r);
    p1->velocity[0] -= dt*diff[0]*mu/m1;
    p1->velocity[1] -= dt*diff[1]*mu/m1;
    p1->velocity[2] -= dt*diff[2]*mu/m1;
    p2->velocity[0] += dt*diff[0]*mu/m2;
    p2->velocity[1] += dt*diff[1]*mu/m2;
    p2->velocity[2] += dt*diff[2]*mu/m2;
}

static void pairwise_gravitational_force_perf_c(const particle *p1, const particle *p2, void *data0, float dt, performance_count *out) {
    (void) p1; (void) p2; (void) data0; (void) dt;
    out->add += 11;
    out->mul += 19;
    out->div += 7;
    out->sqrt += 1;
    out->loads += 15;
    out->stores += 6;
}

static particle_effect_naive pairwise_gravitational_force_effect(float mu) {
    particle_effect_naive result;
    result.particles = 2;
    result.apply.two  = pairwise_gravitational_force_apply;
    result.perf_c.two = pairwise_gravitational_force_perf_c;
    float *data = malloc(1*sizeof(float));
    data[0] = mu;
    result.userdata = data;
    return result;
}

static void pairwise_sphere_collision_apply(particle *p1, particle *p2, void *data0, float dt) {
    (void) dt;
    float *data = data0;
    float radius = data[0];
    float diff[3];
    diff[0] = p2->position[0]-p1->position[0];
    diff[1] = p2->position[1]-p1->position[1];
    diff[2] = p2->position[2]-p1->position[2];
    float r = sqrtf(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
    diff[0] *= 1.0f/r;
    diff[1] *= 1.0f/r;
    diff[2] *= 1.0f/r;
    // project velocities onto diff vector:
    float u1 = p1->velocity[0]*diff[0] + p1->velocity[1]*diff[1] + p1->velocity[2]*diff[2];
    float u2 = p2->velocity[0]*diff[0] + p2->velocity[1]*diff[1] + p2->velocity[2]*diff[2];

    if(r<radius) {
        // within collision range:
        if(u2<u1) {
            float m1 = p1->mass;
            float m2 = p2->mass;
            float restitution = data[1];
            // particles are approaching:
            // calculate new velocities according to collision laws: (in the direction of diff)
            float v1 = (u1*m1+u2*m2-restitution*m2*(u1-u2))/(m1+m2);
            float v2 = (u2*m2+u1*m1-restitution*m1*(u2-u1))/(m1+m2);
            // update velocities:
            p1->velocity[0] += (v1-u1)*diff[0];
            p1->velocity[1] += (v1-u1)*diff[1];
            p1->velocity[2] += (v1-u1)*diff[2];
            p2->velocity[0] += (v2-u2)*diff[0];
            p2->velocity[1] += (v2-u2)*diff[1];
            p2->velocity[2] += (v2-u2)*diff[2];
        }
        // move particles apart, such that they don't intersect (in the direction of diff)
        p1->position[0] -= 0.5f*(radius-r)*diff[0];
        p1->position[1] -= 0.5f*(radius-r)*diff[1];
        p1->position[2] -= 0.5f*(radius-r)*diff[2];
        p2->position[0] += 0.5f*(radius-r)*diff[0];
        p2->position[1] += 0.5f*(radius-r)*diff[1];
        p2->position[2] += 0.5f*(radius-r)*diff[2];
    }
}

static void pairwise_sphere_collision_perf_c(const particle *p1, const particle *p2, void *data0, float dt, performance_count *out) {
    (void) dt;
    float *data = data0;
    float radius = data[0];
    float diff[3];
    diff[0] = p2->position[0]-p1->position[0];
    diff[1] = p2->position[1]-p1->position[1];
    diff[2] = p2->position[2]-p1->position[2];
    float r = sqrtf(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
    diff[0] *= 1.0f/r;
    diff[1] *= 1.0f/r;
    diff[2] *= 1.0f/r;
    // project velocities onto diff vector:
    float u1 = p1->velocity[0]*diff[0] + p1->velocity[1]*diff[1] + p1->velocity[2]*diff[2];
    float u2 = p2->velocity[0]*diff[0] + p2->velocity[1]*diff[1] + p2->velocity[2]*diff[2];
    out->add += 9;
    out->mul += 12;
    out->div += 1;
    out->sqrt += 1;
    out->loads += 13;
    out->cmp += 1;
    if(r<radius) {
        out->cmp += 1;
        if(u2<u1) {
            out->loads += 3;
            out->mul += 12;
            out->add += 15;
            out->stores += 6;
        }
        out->mul += 4;
        out->add += 7;
        out->stores += 6;
    }
}

static particle_effect_naive pairwise_sphere_collision_effect(float radius, float restitution) {
    particle_effect_naive result;
    result.particles = 2;
    result.apply.two  = pairwise_sphere_collision_apply;
    result.perf_c.two = pairwise_sphere_collision_perf_c;
    float *data = malloc(2*sizeof(float));
    data[0] = radius;
    data[1] = restitution;
    result.userdata = data;
    return result;
}

static void newton_step_apply(particle *p, void *data0, float dt) {
    (void) data0;
    p->position[0] += dt*p->velocity[0];
    p->position[1] += dt*p->velocity[1];
    p->position[2] += dt*p->velocity[2];
}

static void newton_step_perf_c(const particle *p, void *data0, float dt, performance_count *out) {
    (void) p; (void) data0; (void) dt;
    out->add += 3;
    out->mul += 3;
    out->loads += 6;
    out->stores += 3;
}

static particle_effect_naive newton_step_effect() {
    particle_effect_naive result;
    result.particles = 1;
    result.apply.one  = newton_step_apply;
    result.perf_c.one = newton_step_perf_c;
    result.userdata = NULL;
    return result;
}

