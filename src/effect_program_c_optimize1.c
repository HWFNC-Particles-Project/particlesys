#include "effect_program_c.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*
 * particle effects consists of a function (apply) and userdata for
 * effect specific data that is passed as second function to apply
 */
typedef struct particle_effect_c_o1_t {
    int particles;
    union {
        void (*one)(particle *, const void *, float);
        void (*two)(particle *, particle *, const void *, float);
    } apply;
    union {
        void (*one)(const particle *, void *, float, performance_count *out);
        void (*two)(const particle *, const particle *, void *, float, performance_count *out);
    } perf_c;
    void *userdata;
} particle_effect_c_o1;


void effect_program_c_o1_compile(effect_program *self, const effect_desc *desc);
void effect_program_c_o1_execute(const effect_program *self, particle_array *arr, float dt);
void effect_program_c_o1_perf_c(const effect_program *self, const particle_array *arr, float dt, performance_count *out);
void effect_program_c_o1_destroy(effect_program *self);

static particle_effect_c_o1 linear_accel_effect(float x, float y, float z);
static particle_effect_c_o1 linear_force_effect(float x, float y, float z);
static particle_effect_c_o1 gravitational_force_effect(float x, float y, float z, float mu);
static particle_effect_c_o1 plane_bounce_effect(float x, float y, float z, float d, float a);
static particle_effect_c_o1 sphere_bounce_effect(float x, float y, float z, float r, float a);
static particle_effect_c_o1 pairwise_gravitational_force_effect(float mu);
static particle_effect_c_o1 pairwise_sphere_collision_effect(float radius, float restitution);
static particle_effect_c_o1 newton_step_effect();

static void linear_accel_apply(particle *p, const void *data0, float dt);
static void linear_force_apply(particle *p, const void *data0, float dt);
static void gravitational_force_apply(particle *p, const void *data0, float dt);
static void plane_bounce_apply(particle *p, const void *data0, float dt);
static void sphere_bounce_apply(particle *p, const void *data0, float dt);
static void newton_step_apply(particle *p, const void *data0, float dt);



int effect_program_create_c_optimze1(effect_program *p) {
    memset(p, 0, sizeof(effect_program));
    p->compile = effect_program_c_o1_compile;
    p->execute = effect_program_c_o1_execute;
    p->perf_c  = effect_program_c_o1_perf_c;
    p->destroy = effect_program_c_o1_destroy;
    p->usr = NULL;
    return 0;
}

void effect_program_c_o1_destroy(effect_program *self) {
    if (self->usr != NULL) {
        // free last compilation result:
        particle_effect_c_o1 *effects = (particle_effect_c_o1 *)self->usr;
        for(size_t j = 0; effects[j].apply.one != NULL; ++j) {
            if (effects[j].userdata != NULL) {
                free(effects[j].userdata);
            }
        }
        free(self->usr);
        self->usr = NULL;
    }
}

void effect_program_c_o1_compile(effect_program *self, const effect_desc *desc) {
    if (self->usr != NULL) {
        effect_program_c_o1_destroy(self);
    }
    // compile ...
    size_t count = desc->size;
    self->usr = malloc((count + 1) * sizeof(particle_effect_c_o1));
    particle_effect_c_o1 *effects = (particle_effect_c_o1 *)self->usr;
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

void effect_program_c_o1_execute(const effect_program *self, particle_array *arr, float dt) {
    if (self->usr == NULL) {
        // fail, no compilation result.
        return;
    }
    particle_effect_c_o1 *effects = (particle_effect_c_o1 *)self->usr;
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

void effect_program_c_o1_perf_c(const effect_program *self, const particle_array *arr, float dt, performance_count *out) {
    // set all counters to zero:
    memset(out, 0, sizeof(*out));
    if (self->usr == NULL) {
        // fail, no compilation result.
        return;
    }
    particle_effect_c_o1 *effects = (particle_effect_c_o1 *)self->usr;
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

static void linear_accel_apply(particle *p, const void *data0, float dt) {
    const float *f_data = (const float *)data0;
    float fd0 = f_data[0];
    float fd1 = f_data[1];
    float fd2 = f_data[2];
    float pv0 = p->velocity[0];
    float pv1 = p->velocity[1];
    float pv2 = p->velocity[2];
    float fd0_dt = dt * fd0;
    float fd1_dt = dt * fd1;
    float fd2_dt = dt * fd2;
    p->velocity[0] = pv0 + fd0_dt;
    p->velocity[1] = pv1 + fd1_dt;
    p->velocity[2] = pv2 + fd2_dt;
}

static void linear_accel_perf_c(const particle *p, void *data0, float dt, performance_count *out) {
    (void) p; (void) data0; (void) dt;
    out->add += 3;
    out->mul += 3;
    out->loads += 6;
    out->stores += 3;
}

static particle_effect_c_o1 linear_accel_effect(float x, float y, float z) {
    particle_effect_c_o1 result;
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

static void linear_force_apply(particle *p, const void *data0, float dt) {
    const float *f_data = (const float *)data0;
    float m = p->mass;
    float k = dt / m;
    float fd0 = f_data[0];
    float fd1 = f_data[1];
    float fd2 = f_data[2];
    float pv0 = p->velocity[0];
    float pv1 = p->velocity[1];
    float pv2 = p->velocity[2];
    float fd0_k = k * fd0;
    float fd1_k = k * fd1;
    float fd2_k = k * fd2;
    p->velocity[0] = pv0 + fd0_k;
    p->velocity[1] = pv1 + fd1_k;
    p->velocity[2] = pv2 + fd2_k;
}

static void linear_force_perf_c(const particle *p, void *data0, float dt, performance_count *out) {
    (void) p; (void) data0; (void) dt;
    out->add += 3;
    out->mul += 3;
    out->div += 1;
    out->loads += 7;
    out->stores += 3;
}

static particle_effect_c_o1 linear_force_effect(float x, float y, float z) {
    particle_effect_c_o1 result;
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

static void gravitational_force_apply(particle *p, const void *data0, float dt) {
    const float *f_data = (const float *)data0;
    float fd0 = f_data[0];
    float fd1 = f_data[1];
    float fd2 = f_data[2];
    float mu = f_data[3];
    float pp0 = p->position[0];
    float pp1 = p->position[1];
    float pp2 = p->position[2];
    float pv0 = p->velocity[0];
    float pv1 = p->velocity[1];
    float pv2 = p->velocity[2];
    float d0 = pp0-fd0;
    float d1 = pp1-fd1;
    float d2 = pp2-fd2;
    float r2 = d0*d0 + d1*d1 + d2*d2;
    float r = sqrt(r2);
    float k_0 = dt * mu;
    float k_1 = r * r2;
    float k_2 = k_0 / k_1;
    float d0_k = k_2 * d0;
    float d1_k = k_2 * d1;
    float d2_k = k_2 * d2;
    p->velocity[0] = pv0 + d0_k;
    p->velocity[1] = pv1 + d1_k;
    p->velocity[2] = pv2 + d2_k;
}

static void gravitational_force_perf_c(const particle *p, void *data0, float dt, performance_count *out) {
    (void) p; (void) data0; (void) dt;
    out->add += 8;
    out->mul += 8;
    out->div += 1;
    out->sqrt += 1;
    out->loads += 10;
    out->stores += 3;
}

static particle_effect_c_o1 gravitational_force_effect(float x, float y, float z, float mu) {
    particle_effect_c_o1 result;
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

static void plane_bounce_apply(particle *p, const void *data0, float dt) {
    (void) dt;
    const float *f_data = (const float *)data0;
    float a = f_data[0];
    float b = f_data[1];
    float c = f_data[2];
    float d = f_data[3];
    float dec = f_data[4];
    float pp0 = p->position[0];
    float pp1 = p->position[1];
    float pp2 = p->position[2];
    float pv0 = p->velocity[0];
    float pv1 = p->velocity[1];
    float pv2 = p->velocity[2];
    float a_pp0 = a * pp0;
    float b_pp1 = b * pp1;
    float c_pp2 = c * pp2;
    float d_a_pp0 = d - a_pp0;
    float b_pp1_c_pp2 = b_pp1 + c_pp2;
    float a_pv0 = a * pv0;
    float b_pv1 = b * pv1;
    float c_pv2 = c * pv2;
    float a_pv0_b_pv1 = a_pv0 +  b_pv1;
    float d_dist = d_a_pp0 - b_pp1_c_pp2;
    float vnormal = a_pv0_b_pv1 + c_pv2;
    //float dist = a_pp0 + b_pp1 + c_pp2;
    if(d_dist > 0.0f && vnormal < 0.0f) {
        // we are behind plane and velocity is away from the back of the plane
        float vnormal_a = vnormal * a;
        float vnormal_b = vnormal * b;
        float vnormal_c = vnormal * c;
        //float d_dist = d - dist;
        float pv0_1 = pv0 - vnormal_a-vnormal_a;
        float pv1_1 = pv1 - vnormal_b-vnormal_b;
        float pv2_1 = pv2 - vnormal_c-vnormal_c;
        float d_dist_a = d_dist * a;
        float d_dist_b = d_dist * b;
        float d_dist_c = d_dist * c;
        p->position[0] = pp0 + d_dist_a;
        p->position[1] = pp1 + d_dist_b;
        p->position[2] = pp2 + d_dist_c;
        p->velocity[0] = pv0_1 * dec;
        p->velocity[1] = pv1_1 * dec;
        p->velocity[2] = pv2_1 * dec;
    }
}

static void plane_bounce_perf_c(const particle *p, void *data0, float dt, performance_count *out) {
    (void) dt;
    float *data = data0;
    float dist =    data[0]*p->position[0] + data[1]*p->position[1] + data[2]*p->position[2];
    float vnormal = data[0]*p->velocity[0] + data[1]*p->velocity[1] + data[2]*p->velocity[2];
    float d = data[3];
    out->add += 5;
    out->mul += 6;
    out->loads += 11;
    out->cmp += 1;
    if(dist<d) {
        out->cmp += 1;
        if (vnormal<0.0f) {
            // branch is taken:
            out->add += 9;
            out->mul += 9;
            out->loads += 0;
            out->stores += 6;
        }
    } else {
        // branch not taken.
    }
}

static particle_effect_c_o1 plane_bounce_effect(float x, float y, float z, float d, float a) {
    particle_effect_c_o1 result;
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

static void sphere_bounce_apply(particle *p, const void *data0, float dt) {
    (void) dt;
    const float *f_data = (const float *)data0;
    float fd0 = f_data[0];
    float fd1 = f_data[1];
    float fd2 = f_data[2];
    float d   = f_data[3];
    float dec = f_data[4];
    float pp0 = p->position[0];
    float pp1 = p->position[1];
    float pp2 = p->position[2];
    float pv0 = p->velocity[0];
    float pv1 = p->velocity[1];
    float pv2 = p->velocity[2];
    float d2 = d * d;
    float n0  = pp0-fd0;
    float n1  = pp1-fd1;
    float n2  = pp2-fd2;
    float r2 = n0*n0 + n1*n1 + n2*n2;
    float vnormal = n0*pv0 + n1*pv1 + n2*pv2;
    float r2_reci = 1.0f / r2;
    if(d2 > r2 && vnormal < 0.0f) {
        // inside sphere and going further inside sphere:
        float r_reci = sqrtf(r2_reci);
        
        float vnormal2_r2 = 2.0f * vnormal * r2_reci;
        
        float pv0_1 = pv0 - vnormal2_r2*n0;
        float pv1_1 = pv1 - vnormal2_r2*n1;
        float pv2_1 = pv2 - vnormal2_r2*n2;
        
        float d_r_r_reci = d * r_reci - 1.0f;

        p->velocity[0] = pv0_1 * dec;
        p->velocity[1] = pv1_1 * dec;
        p->velocity[2] = pv2_1 * dec;

        p->position[0] = pp0 + d_r_r_reci*n0;
        p->position[1] = pp1 + d_r_r_reci*n1;
        p->position[2] = pp2 + d_r_r_reci*n2;
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
    out->mul += 7;
    out->div += 1;
    out->loads += 11;
    out->cmp += 1;
    if(r<d) {
        out->cmp += 1;
        if (vnormal<0.0f) {
            // branch is taken:
            out->add += 7;
            out->sqrt += 1;
            out->mul += 12;
            out->stores += 6;
        }
    } else {
        // branch not taken.
    }
}

static particle_effect_c_o1 sphere_bounce_effect(float x, float y, float z, float r, float a) {
    particle_effect_c_o1 result;
    result.particles = 1;
    result.apply.one =  sphere_bounce_apply;
    result.perf_c.one = sphere_bounce_perf_c;
    float *data = malloc(5*sizeof(float));
    data[0] = x;
    data[1] = y;
    data[2] = z;
    data[3] = r;
    //data[4] = r * r;
    data[4] = a;
    result.userdata = data;
    return result;
}

static void pairwise_gravitational_force_apply(particle *p1, particle *p2, const void *data0, float dt) {
    const float *f_data = (const float *)data0;
    float mu = f_data[0];
    float pp10 = p1->position[0];
    float pp11 = p1->position[1];
    float pp12 = p1->position[2];
    float pp20 = p2->position[0];
    float pp21 = p2->position[1];
    float pp22 = p2->position[2];
    float pv10 = p1->velocity[0];
    float pv11 = p1->velocity[1];
    float pv12 = p1->velocity[2];
    float pv20 = p2->velocity[0];
    float pv21 = p2->velocity[1];
    float pv22 = p2->velocity[2];
    float m1 = p1->mass;
    float m2 = p2->mass;
    float dt_mu = dt * mu;
    float d0 = pp20-pp10;
    float d1 = pp21-pp11;
    float d2 = pp22-pp12;
    float r2 = d0*d0 + d1*d1 + d2*d2;
    float r4 = r2 * r2;
    float r = sqrt(r2);
    // is r^3 because we need to normalize diff.
    float mu_0 = dt_mu / r4;
    float mu_1 = mu_0 * r;
    float mu_m2 = mu_1 * m2;
    float mu_m1 = mu_1 * m1;
    float d0_1 = d0 * mu_m2;
    float d1_1 = d1 * mu_m2;
    float d2_1 = d2 * mu_m2;
    float d0_2 = d0 * mu_m1;
    float d1_2 = d1 * mu_m1;
    float d2_2 = d2 * mu_m1;
    p1->velocity[0] = pv10 - d0_1;
    p1->velocity[1] = pv11 - d1_1;
    p1->velocity[2] = pv12 - d2_1;
    p2->velocity[0] = pv20 + d0_2;
    p2->velocity[1] = pv21 + d1_2;
    p2->velocity[2] = pv22 + d2_2;
}

static void pairwise_gravitational_force_perf_c(const particle *p1, const particle *p2, void *data0, float dt, performance_count *out) {
    (void) p1; (void) p2; (void) data0; (void) dt;
    out->add += 11;
    out->mul += 14;
    out->div += 1;
    out->sqrt += 1;
    out->loads += 15;
    out->stores += 6;
}

static particle_effect_c_o1 pairwise_gravitational_force_effect(float mu) {
    particle_effect_c_o1 result;
    result.particles = 2;
    result.apply.two  = pairwise_gravitational_force_apply;
    result.perf_c.two = pairwise_gravitational_force_perf_c;
    float *data = malloc(1*sizeof(float));
    data[0] = mu;
    result.userdata = data;
    return result;
}

static void pairwise_sphere_collision_apply(particle *p1, particle *p2, const void *data0, float dt) {
    (void) dt;
    const float *f_data = (const float *)data0;
    float radius = f_data[0];
    float pp10 = p1->position[0];
    float pp11 = p1->position[1];
    float pp12 = p1->position[2];
    float pp20 = p2->position[0];
    float pp21 = p2->position[1];
    float pp22 = p2->position[2];
    float pv10 = p1->velocity[0];
    float pv11 = p1->velocity[1];
    float pv12 = p1->velocity[2];
    float pv20 = p2->velocity[0];
    float pv21 = p2->velocity[1];
    float pv22 = p2->velocity[2];
    float d0 = pp20 - pp10;
    float d1 = pp21 - pp11;
    float d2 = pp22 - pp12;
    float radius2 = radius * radius;
    float r2 = d0 * d0 + d1 * d1 + d2 * d2;
    // project velocities onto diff vector:
    float pv10d0 = pv10*d0;
    float pv11d1 = pv11*d1;
    float pv12d2 = pv12*d2;
    float pv20d0 = pv20*d0;
    float pv21d1 = pv21*d1;
    float pv22d2 = pv22*d2;
    float u1_0 = pv10d0 + pv11d1 + pv12d2;
    float u2_0 = pv20d0 + pv21d1 + pv22d2;
    float u1_u2_0 = (pv10d0-pv20d0) + (pv11d1-pv21d1) + (pv12d2-pv22d2);
    //float u1_u2_0 = u1_0 - u2_0;
    if(r2 < radius2) {
        // within collision range:
        float r = sqrt(r2);
        float r_reci = 1.0f / r;
        float radius_r = radius - r;
        float radius_r05 = 0.5f * radius_r;
        float d0_1 = d0 * r_reci;
        float d1_1 = d1 * r_reci;
        float d2_1 = d2 * r_reci;
        if(u1_u2_0 > 0.0f) {
            float m1 = p1->mass;
            float m2 = p2->mass;
            float restitution = f_data[1];
            float m12_reci = 1.0f/(m1+m2);
            float u1 = u1_0 * r_reci;
            float u2 = u2_0 * r_reci;
            float um1pum2 = u1*m1+u2*m2;
            float u1_u2 = u1_u2_0 * r_reci;
            // particles are approaching:
            // calculate new velocities according to collision laws: (in the direction of diff)
            float v1 = (um1pum2-restitution*m2*u1_u2)*m12_reci;
            float v2 = (um1pum2+restitution*m1*u1_u2)*m12_reci;
            // update velocities:
            float v1_u1 = v1 - u1;
            float v2_u2 = v2 - u2;
            float d0_1_v1_u1 = d0_1 * v1_u1;
            float d1_1_v1_u1 = d1_1 * v1_u1;
            float d2_1_v1_u1 = d2_1 * v1_u1;
            float d0_1_v2_u2 = d0_1 * v2_u2;
            float d1_1_v2_u2 = d1_1 * v2_u2;
            float d2_1_v2_u2 = d2_1 * v2_u2;
            p1->velocity[0] = pv10 + d0_1_v1_u1;
            p1->velocity[1] = pv11 + d1_1_v1_u1;
            p1->velocity[2] = pv12 + d2_1_v1_u1;
            p2->velocity[0] = pv20 + d0_1_v2_u2;
            p2->velocity[1] = pv21 + d1_1_v2_u2;
            p2->velocity[2] = pv22 + d2_1_v2_u2;
        }
        // move particles apart, such that they don't intersect (in the direction of diff)
        float d0_1_rr05 = d0_1*radius_r05;
        float d1_1_rr05 = d1_1*radius_r05;
        float d2_1_rr05 = d2_1*radius_r05;
        p1->position[0] = pp10 - d0_1_rr05;
        p1->position[1] = pp11 - d1_1_rr05;
        p1->position[2] = pp12 - d2_1_rr05;
        p2->position[0] = pp20 + d0_1_rr05;
        p2->position[1] = pp21 + d1_1_rr05;
        p2->position[2] = pp22 + d2_1_rr05;
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
    float r = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
    diff[0] *= 1.0f/r;
    diff[1] *= 1.0f/r;
    diff[2] *= 1.0f/r;
    // project velocities onto diff vector:
    float u1 = p1->velocity[0]*diff[0] + p1->velocity[1]*diff[1] + p1->velocity[2]*diff[2];
    float u2 = p2->velocity[0]*diff[0] + p2->velocity[1]*diff[1] + p2->velocity[2]*diff[2];
    out->add += 14;
    out->mul += 10;
    out->loads += 13;
    out->cmp += 1;
    if(r<radius) {
        out->div += 1;
        out->sqrt += 1;
        out->add += 1;
        out->mul += 4;
        out->cmp += 1;
        if(u2<u1) {
            out->loads += 3;
            out->div += 1;
            out->mul += 17;
            out->add += 12;
            out->stores += 6;
        }
        out->mul += 3;
        out->add += 6;
        out->stores += 6;
    }
}

static particle_effect_c_o1 pairwise_sphere_collision_effect(float radius, float restitution) {
    particle_effect_c_o1 result;
    result.particles = 2;
    result.apply.two  = pairwise_sphere_collision_apply;
    result.perf_c.two = pairwise_sphere_collision_perf_c;
    float *data = malloc(2*sizeof(float));
    data[0] = radius;
    data[1] = restitution;
    result.userdata = data;
    return result;
}

static void newton_step_apply(particle *p, const void *data0, float dt) {
    (void) data0;
    float pp0 = p->position[0];
    float pp1 = p->position[1];
    float pp2 = p->position[2];
    float pv0 = p->velocity[0];
    float pv1 = p->velocity[1];
    float pv2 = p->velocity[2];
    float pv0_1 = pv0 * dt;
    float pv1_1 = pv1 * dt;
    float pv2_1 = pv2 * dt;
    p->position[0] = pp0 + pv0_1;
    p->position[1] = pp1 + pv1_1;
    p->position[2] = pp2 + pv2_1;
}

static void newton_step_perf_c(const particle *p, void *data0, float dt, performance_count *out) {
    (void) p; (void) data0; (void) dt;
    out->add += 3;
    out->mul += 3;
    out->loads += 6;
    out->stores += 3;
}

static particle_effect_c_o1 newton_step_effect() {
    particle_effect_c_o1 result;
    result.particles = 1;
    result.apply.one  = newton_step_apply;
    result.perf_c.one = newton_step_perf_c;
    result.userdata = NULL;
    return result;
}

