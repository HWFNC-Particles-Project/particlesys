#include "effect_program_c.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <xmmintrin.h>
#include <smmintrin.h>
#include "malloc_align.h"

/*
 * particle effects consists of a function (apply) and userdata for
 * effect specific data that is passed as second function to apply
 */
typedef struct particle_effect_c_o2_t {
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
    void *p_to_free;
} particle_effect_c_o2;


void effect_program_c_o2_compile(effect_program *self, const effect_desc *desc);
void effect_program_c_o2_execute(const effect_program *self, particle_array *arr, float dt);
void effect_program_c_o2_perf_c(const effect_program *self, const particle_array *arr, float dt, performance_count *out);
void effect_program_c_o2_destroy(effect_program *self);

static particle_effect_c_o2 linear_accel_effect(float x, float y, float z);
static particle_effect_c_o2 linear_force_effect(float x, float y, float z);
static particle_effect_c_o2 gravitational_force_effect(float x, float y, float z, float mu);
static particle_effect_c_o2 plane_bounce_effect(float x, float y, float z, float d, float a);
static particle_effect_c_o2 sphere_bounce_effect(float x, float y, float z, float r, float a);
static particle_effect_c_o2 pairwise_gravitational_force_effect(float mu);
static particle_effect_c_o2 pairwise_sphere_collision_effect(float radius, float restitution);
static particle_effect_c_o2 newton_step_effect();

static void linear_accel_apply(particle *p, const void *data0, float dt);
static void linear_force_apply(particle *p, const void *data0, float dt);
static void gravitational_force_apply(particle *p, const void *data0, float dt);
static void plane_bounce_apply(particle *p, const void *data0, float dt);
static void sphere_bounce_apply(particle *p, const void *data0, float dt);
static void newton_step_apply(particle *p, const void *data0, float dt);



int effect_program_create_c_optimze2(effect_program *p) {
    memset(p, 0, sizeof(effect_program));
    p->compile = effect_program_c_o2_compile;
    p->execute = effect_program_c_o2_execute;
    p->perf_c  = effect_program_c_o2_perf_c;
    p->destroy = effect_program_c_o2_destroy;
    p->usr = NULL;
    return 0;
}

void effect_program_c_o2_destroy(effect_program *self) {
    if (self->usr != NULL) {
        // free last compilation result:
        particle_effect_c_o2 *effects = (particle_effect_c_o2 *)self->usr;
        for(size_t j = 0; effects[j].apply.one != NULL; ++j) {
            if (effects[j].p_to_free != NULL) {
                free(effects[j].p_to_free);
            }
        }
        free(self->usr);
        self->usr = NULL;
    }
}

void effect_program_c_o2_compile(effect_program *self, const effect_desc *desc) {
    if (self->usr != NULL) {
        effect_program_c_o2_destroy(self);
    }
    // compile ...
    size_t count = desc->size;
    self->usr = malloc((count + 1) * sizeof(particle_effect_c_o2));
    particle_effect_c_o2 *effects = (particle_effect_c_o2 *)self->usr;
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

void effect_program_c_o2_execute(const effect_program *self, particle_array *arr, float dt) {
    if (self->usr == NULL) {
        // fail, no compilation result.
        return;
    }
    particle_effect_c_o2 *effects = (particle_effect_c_o2 *)self->usr;
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

void effect_program_c_o2_perf_c(const effect_program *self, const particle_array *arr, float dt, performance_count *out) {
    // set all counters to zero:
    memset(out, 0, sizeof(*out));
    if (self->usr == NULL) {
        // fail, no compilation result.
        return;
    }
    particle_effect_c_o2 *effects = (particle_effect_c_o2 *)self->usr;
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
    __m128 dtv  = _mm_set_ps1(dt);
    __m128 a    = _mm_load_ps(&f_data[0]);
    __m128 v    = _mm_load_ps(&p->velocity[0]);
    __m128 dt_a = _mm_mul_ps(dtv, a);
    __m128 v_1  = _mm_add_ps(v, dt_a);
    _mm_store_ps(&p->velocity[0], v_1);
}

static void linear_accel_perf_c(const particle *p, void *data0, float dt, performance_count *out) {
    (void) p; (void) data0; (void) dt;
    out->add += 4;
    out->mul += 4;
    out->loads += 8;
    out->stores += 4;
}

static particle_effect_c_o2 linear_accel_effect(float x, float y, float z) {
    particle_effect_c_o2 result;
    result.particles = 1;
    result.apply.one  = linear_accel_apply;
    result.perf_c.one = linear_accel_perf_c;
    float *data = malloc_align(4*sizeof(float), 16, &result.p_to_free);
    data[0] = x;
    data[1] = y;
    data[2] = z;
    data[3] = 0;
    result.userdata = data;
    return result;
}

static void linear_force_apply(particle *p, const void *data0, float dt) {
    const float *f_data = (const float *)data0;
    __m128 dtv  = _mm_set_ps1(dt);
    __m128 m    = _mm_load_ps1(&p->mass);
    __m128 f    = _mm_load_ps(&f_data[0]);
    __m128 v    = _mm_load_ps(&p->velocity[0]);
    __m128 dt_m = _mm_div_ps(dtv, m);
    __m128 vd   = _mm_mul_ps(dt_m, f);
    __m128 v_1  = _mm_add_ps(v, vd);
    _mm_store_ps(&p->velocity[0], v_1);
}

static void linear_force_perf_c(const particle *p, void *data0, float dt, performance_count *out) {
    (void) p; (void) data0; (void) dt;
    out->add += 4;
    out->mul += 4;
    out->div += 4;
    out->loads += 9;
    out->stores += 4;
}

static particle_effect_c_o2 linear_force_effect(float x, float y, float z) {
    particle_effect_c_o2 result;
    result.particles = 1;
    result.apply.one  = linear_force_apply;
    result.perf_c.one = linear_force_perf_c;
    float *data = malloc_align(4*sizeof(float), 16, &result.p_to_free);
    data[0] = x;
    data[1] = y;
    data[2] = z;
    data[3] = 0;
    result.userdata = data;
    return result;
}

static void print__m128(__m128 d) {
    float data[4];
    _mm_storeu_ps(data, d);
    printf("%12.4e, %12.4e, %12.4e, %12.4e\n", data[0], data[1], data[2], data[3]);
}

static void gravitational_force_apply(particle *p, const void *data0, float dt) {
    const float *f_data = (const float *)data0;
    __m128 dtv  = _mm_set_ps1(dt);
    __m128 cntr = _mm_load_ps(&f_data[0]);
    __m128 mu   = _mm_load_ps(&f_data[4]);
    __m128 pp   = _mm_load_ps(&p->position[0]);
    __m128 v    = _mm_load_ps(&p->velocity[0]);
    //__m128 s15  = _mm_set1_ps(1.5f);
    //__m128 s05  = _mm_set1_ps(-0.5f);
    //__m128 z    = _mm_setzero_ps();      // 1 cycle
    __m128 d    = _mm_sub_ps(pp, cntr);  // 3 cycles    k       -> d = pp - cntr
    __m128 d2   = _mm_mul_ps(d, d);      // 5 cycles    k       -> d2 = d * d
    __m128 mu_dt= _mm_mul_ps(mu, dtv);   // 5 cycles
    __m128 d2_3 = _mm_movehl_ps(d2, d2); // 1 cycle     k       -> d2_3 = [d2[2]       xx    d2[2] xx   ]
    __m128 d2_2 = _mm_shuffle_ps(d2, d2, 1);// 1 cycle          -> d2_2 = [d2[1]       d2[0] d2[0] d2[0]]
    __m128 d2_1 = _mm_add_ss(d2,   d2_3);// 3 cycles    k       -> d2_1 = [d2[0]+d2[2] d2[1] d2[2] xx   ]
    __m128 r2   = _mm_add_ss(d2_2, d2_1);// 3 cycles    k       -> r2   = [r2          d2[0] d2[0] d2[0]]
    __m128 d_mu_dt= _mm_mul_ps(d, mu_dt);// 5 cycles
    __m128 r2_a = _mm_shuffle_ps(r2, r2, 0);// 1 cycle          -> r2_a = [r2 r2 r2 r2]
    
    __m128 r    = _mm_sqrt_ps(r2_a);     // 20 cycles   k
    __m128 r4   = _mm_mul_ps(r2_a, r2_a);// 5 cycles    i
    // reciprocal with one newton step for accuracy:
    __m128 r4_ra= _mm_rcp_ps(r4);        // 3 cycles    i
    __m128 r4_r2= _mm_add_ps(r4_ra, r4_ra);// 3 cycles
    __m128 r8_ra= _mm_mul_ps(r4_ra, r4_ra);// 5 cycles  i
    __m128 r4_r0= _mm_mul_ps(r8_ra, r4);   // 5 cycles  i
    __m128 r4_r = _mm_sub_ps(r4_r2, r4_r0);// 3 cycles  i
    
    __m128 d_r  = _mm_mul_ps(d_mu_dt, r);  // 5 cycles    k
    __m128 d_k_2= _mm_mul_ps(d_r, r4_r);  // 5 cycles    k
    __m128 v_1  = _mm_add_ps(v, d_k_2);   // 3 cycles    k
    _mm_store_ps(&p->velocity[0], v_1);
}
    /* base case, but sqrt and div do not run in parallel.
    __m128 r    = _mm_sqrt_ps(r2_a);     // 20 cycles   k
    __m128 r4   = _mm_mul_ps(r2_a, r2_a);// 5 cycles    i
    __m128 k_1  = _mm_div_ps(mu_dt, r4); // 14 cycles   i
    __m128 k_2  = _mm_mul_ps(k_1, r);    // 5 cycles    k
    */
    /* very poor accuracy but fast:
    __m128 r_rec= _mm_rsqrt_ps(r2_a);       // 5 cycles    k
    __m128 k_0  = _mm_mul_ps(mu_dt, r_rec); // 5 cycles    k
    __m128 k_1  = _mm_mul_ps(r_rec, r_rec); // 5 cycles
    __m128 k_2  = _mm_mul_ps(k_0, k_1);     // 5 cycles    k+1
    */
    /* not faster but accurate
    __m128 r4   = _mm_mul_ps(r2_a, r2_a);// 5 cycles  
    __m128 s05r2= _mm_mul_ps(r2_a, s05);
    __m128 r_r_0= _mm_rsqrt_ps(r2_a);    // 5 cycles   k
    __m128 r_r_2= _mm_mul_ps(r_r_0, r_r_0);// 5 cycles k
    __m128 s15rr= _mm_mul_ps(s15, r_r_0);
    __m128 s05r0= _mm_mul_ps(s05r2, r_r_0);
    __m128 s05_0= _mm_mul_ps(s05r0, r_r_2);// 5 cycles k
    __m128 r_r  = _mm_add_ps(s15rr, s05_0);// 3 cycles k
    __m128 r    = _mm_mul_ps(r_r, r2_a); // 5 cycles   
    __m128 k_1  = _mm_div_ps(mu_dt, r4); // 14 cycles  k
    __m128 k_2  = _mm_mul_ps(k_1, r);    // 5 cycles   k
    */
    //__m128 d_k_2= _mm_mul_ps(d, k_2);    // 5 cycles    k
    //__m128 v_1  = _mm_add_ps(v, d_k_2);  // 3 cycles    k

static void gravitational_force_perf_c(const particle *p, void *data0, float dt, performance_count *out) {
    (void) p; (void) data0; (void) dt;
    out->add += 18;
    out->mul += 32;
    out->div += 4;
    out->sqrt += 4;
    out->loads += 16;
    out->stores += 4;
}

static particle_effect_c_o2 gravitational_force_effect(float x, float y, float z, float mu) {
    particle_effect_c_o2 result;
    result.particles = 1;
    result.apply.one  = gravitational_force_apply;
    result.perf_c.one = gravitational_force_perf_c;
    float *data = malloc_align(8*sizeof(float), 16, &result.p_to_free);
    data[0] = x;
    data[1] = y;
    data[2] = z;
    data[3] = 0;
    data[4] = mu;
    data[5] = mu;
    data[6] = mu;
    data[7] = 0;
    result.userdata = data;
    result.p_to_free = data;
    return result;
}

static void plane_bounce_apply(particle *p, const void *data0, float dt) {
    (void) dt;
    const float *f_data = (const float *)data0;
    __m128 z    = _mm_setzero_ps();
    __m128 nv   = _mm_load_ps(&f_data[0]);
    __m128 pp   = _mm_load_ps(&p->position[0]);
    __m128 v    = _mm_load_ps(&p->velocity[0]);
    __m128 pp1  = _mm_blend_ps(pp, _mm_set_ps1(-1.0f), 0b1000);
    // multiply vectors:
    __m128 pd   = _mm_mul_ps(pp1, nv);
    __m128 vd   = _mm_mul_ps(v  , nv);
    // horizontal add them:
    __m128 vd_3 = _mm_movehl_ps(vd, vd); // 1 cycle     k       -> vd_3 = [vd[2]       xx          vd[2]       xx   ]
    __m128 pd34 = _mm_movehl_ps(pd, pd); // 1 cycle     k       -> pd34 = [pd[2]       pd[3]       pd[2]       pd[3]]
    __m128 vd_2 = _mm_shuffle_ps(vd, vd, 1);// 1 cycle          -> vd_2 = [vd[1]       vd[0]       vd[0]       vd[0]]
    __m128 pd_1 = _mm_add_ps(pd,   pd34);// 3 cycles    k       -> pd_1 = [pd[0]+pd[2] pd[1]+pd[3] pd[2]       pd[3]]
    __m128 vd_1 = _mm_add_ss(vd,   vd_3);// 3 cycles    k       -> vd_1 = [vd[0]+vd[2] vd[1]       vd[2]       xx   ]
    __m128 pd_2 = _mm_shuffle_ps(pd_1, pd_1, 1);// 1 cycle      -> pd_2 = [pd[1]+pd[3] pd[0]+pd[2] pd[0]+pd[2] pd[0]+pd[2]]
    __m128 ddst = _mm_add_ps(pd_2, pd_1);// 3 cycles    k       -> vnrm = [ddst        ddst        xx          xx   ]
    __m128 vnrm = _mm_add_ss(vd_2, vd_1);// 3 cycles    k       -> vnrm = [vnrm        vd[0]       vd[0]       vd[0]]
    
    /*union {
        __m128 i;
        uint32_t f[4];
    } brnch = {.i = _mm_and_ps(_mm_cmpgt_ss(z, ddst), _mm_cmpgt_ss(z, vnrm))};
    
    if(brnch.f[0]) {*/
    float d_dist; float vnormal;
    _mm_store_ss(&d_dist, ddst); _mm_store_ss(&vnormal, vnrm);
    if (d_dist < 0.0f && vnormal < 0.0f) {
        // we are behind plane and velocity is away from the back of the plane
        __m128 nv_z   = _mm_blend_ps(nv, z, 0b1000);
        __m128 ddst_a = _mm_movelh_ps(ddst, ddst);    // 1 cycle          -> ddst_a = [ddst ddst ddst ddst]
        __m128 vnrm_a = _mm_shuffle_ps(vnrm, vnrm, 0);// 1 cycle          -> vnrm_a = [vnrm vnrm vnrm vnrm]
        __m128 dec    = _mm_load_ps(&f_data[4]);
        __m128 vnrm_n = _mm_mul_ps(vnrm_a, nv_z);
        __m128 ddst_n = _mm_mul_ps(ddst_a, nv_z);
        __m128 v_1    = _mm_sub_ps(v, vnrm_n);
        __m128 v_2    = _mm_sub_ps(v_1, vnrm_n);
        __m128 pp_1   = _mm_sub_ps(pp, ddst_n);
        __m128 v_dec  = _mm_mul_ps(v_2, dec);
        _mm_store_ps(&p->position[0], pp_1);
        _mm_store_ps(&p->velocity[0], v_dec);
    }
}

static void plane_bounce_perf_c(const particle *p, void *data0, float dt, performance_count *out) {
    (void) dt;
    float *data = data0;
    float dist =    data[0]*p->position[0] + data[1]*p->position[1] + data[2]*p->position[2];
    float vnormal = data[0]*p->velocity[0] + data[1]*p->velocity[1] + data[2]*p->velocity[2];
    float d = data[3];
    out->add += 10;
    out->mul += 8;
    out->loads += 12;
    out->cmp += 1;
    if(dist<d) {
        out->cmp += 1;
        if (vnormal<0.0f) {
            // branch is taken:
            out->add += 12;
            out->mul += 12;
            out->loads += 4;
            out->stores += 8;
        }
    } else {
        // branch not taken.
    }
}

static particle_effect_c_o2 plane_bounce_effect(float x, float y, float z, float d, float a) {
    particle_effect_c_o2 result;
    result.particles = 1;
    result.apply.one =  plane_bounce_apply;
    result.perf_c.one = plane_bounce_perf_c;
    float *data = malloc_align(8*sizeof(float), 16, &result.p_to_free);
    float r = sqrtf(x*x + y*y + z*z);
    data[0] = x/r;
    data[1] = y/r;
    data[2] = z/r;
    data[3] = d/r;
    data[4] = a;
    data[5] = a;
    data[6] = a;
    data[7] = 1;
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

static particle_effect_c_o2 sphere_bounce_effect(float x, float y, float z, float r, float a) {
    particle_effect_c_o2 result;
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
    result.p_to_free = data;
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
    float r = sqrtf(r2);
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

static particle_effect_c_o2 pairwise_gravitational_force_effect(float mu) {
    particle_effect_c_o2 result;
    result.particles = 2;
    result.apply.two  = pairwise_gravitational_force_apply;
    result.perf_c.two = pairwise_gravitational_force_perf_c;
    float *data = malloc(1*sizeof(float));
    data[0] = mu;
    result.userdata = data;
    result.p_to_free = data;
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
        float r = sqrtf(r2);
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
    float r = sqrtf(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
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

static particle_effect_c_o2 pairwise_sphere_collision_effect(float radius, float restitution) {
    particle_effect_c_o2 result;
    result.particles = 2;
    result.apply.two  = pairwise_sphere_collision_apply;
    result.perf_c.two = pairwise_sphere_collision_perf_c;
    float *data = malloc(2*sizeof(float));
    data[0] = radius;
    data[1] = restitution;
    result.userdata = data;
    result.p_to_free = data;
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

static particle_effect_c_o2 newton_step_effect() {
    particle_effect_c_o2 result;
    result.particles = 1;
    result.apply.one  = newton_step_apply;
    result.perf_c.one = newton_step_perf_c;
    result.userdata = NULL;
    result.p_to_free = NULL;
    return result;
}
