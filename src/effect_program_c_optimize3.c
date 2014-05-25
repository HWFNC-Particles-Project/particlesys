#include "effect_program_c.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <xmmintrin.h>
#include <smmintrin.h>
#include "malloc_align.h"

#define M 8

typedef struct reordered_particles_t {
    float position[3][4];
    float     mass[1][4];
    float velocity[3][4];
    float   charge[1][4];
} reordered_particles;

/*
 * particle effects consists of a function (apply) and userdata for
 * effect specific data that is passed as second function to apply
 */
typedef struct particle_effect_c_o3_t {
    int particles;
    union {
        void (*one)(reordered_particles *, const void *, float);
        void (*two)(reordered_particles *, reordered_particles *, const void *, float);
    } apply;
    union {
        void (*one)(const reordered_particles *, void *, float, performance_count *out);
        void (*two)(const reordered_particles *, const reordered_particles *, void *, float, performance_count *out);
    } perf_c;
    void *userdata;
    void *p_to_free;
} particle_effect_c_o3;


void effect_program_c_o3_compile(effect_program *self, const effect_desc *desc);
void effect_program_c_o3_execute(const effect_program *self, particle_array *arr, float dt);
void effect_program_c_o3_perf_c(const effect_program *self, const particle_array *arr, float dt, performance_count *out);
void effect_program_c_o3_destroy(effect_program *self);

static particle_effect_c_o3 linear_accel_effect(float x, float y, float z);
static particle_effect_c_o3 linear_force_effect(float x, float y, float z);
static particle_effect_c_o3 gravitational_force_effect(float x, float y, float z, float mu);
static particle_effect_c_o3 plane_bounce_effect(float x, float y, float z, float d, float a);
static particle_effect_c_o3 sphere_bounce_effect(float x, float y, float z, float r, float a);
static particle_effect_c_o3 pairwise_gravitational_force_effect(float mu);
static particle_effect_c_o3 pairwise_sphere_collision_effect(float radius, float restitution);
static particle_effect_c_o3 newton_step_effect();


int effect_program_create_c_optimze3(effect_program *p) {
    memset(p, 0, sizeof(effect_program));
    p->compile = effect_program_c_o3_compile;
    p->execute = effect_program_c_o3_execute;
    p->perf_c  = effect_program_c_o3_perf_c;
    p->destroy = effect_program_c_o3_destroy;
    p->usr = NULL;
    return 0;
}

void effect_program_c_o3_destroy(effect_program *self) {
    if (self->usr != NULL) {
        // free last compilation result:
        particle_effect_c_o3 *effects = (particle_effect_c_o3 *)self->usr;
        for(size_t j = 0; effects[j].apply.one != NULL; ++j) {
            if (effects[j].p_to_free != NULL) {
                free(effects[j].p_to_free);
            }
        }
        free(self->usr);
        self->usr = NULL;
    }
}

void effect_program_c_o3_compile(effect_program *self, const effect_desc *desc) {
    if (self->usr != NULL) {
        effect_program_c_o3_destroy(self);
    }
    // compile ...
    size_t count = desc->size;
    self->usr = malloc((count + 1) * sizeof(particle_effect_c_o3));
    particle_effect_c_o3 *effects = (particle_effect_c_o3 *)self->usr;
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
                printf("error: sphere bounce not implemented!\n");
                /*effects[i] = sphere_bounce_effect(      el->float_usr[0],
                                                        el->float_usr[1],
                                                        el->float_usr[2],
                                                        el->float_usr[3],
                                                        el->float_usr[4]);*/
                effects[i] = linear_accel_effect(0,0,0);
                break;
            case EFFECT_TYPE_NEWTON_STEP:
                effects[i] = newton_step_effect();
                break;

            case EFFECT_TYPE_GRAVITY_FORCE:
                printf("error: gravity force not implemented!\n");
                //effects[i] = pairwise_gravitational_force_effect(el->float_usr[0]);
                effects[i] = linear_accel_effect(0,0,0);
                break;

            case EFFECT_TYPE_SPHERE_COLLISION:
                printf("error: sphere collision not implemented!\n");
                //effects[i] = pairwise_sphere_collision_effect(  el->float_usr[0],
                //                                                el->float_usr[1]);
                effects[i] = linear_accel_effect(0,0,0);
                break;
            default:
                break; // unhandled effect
        }
    }
    effects[count].apply.one = NULL;
    effects[count].userdata = NULL;
}

void effect_program_c_o3_execute(const effect_program *self, particle_array *arr, float dt) {
    if (self->usr == NULL) {
        // fail, no compilation result.
        return;
    }
    particle_effect_c_o3 *effects = (particle_effect_c_o3 *)self->usr;
    // create working array:
    void *warray1_mem = NULL;
    void *warray2_mem = NULL;
    reordered_particles *warray1 = malloc_align(M * sizeof(reordered_particles), 16, &warray1_mem);
    reordered_particles *warray2 = malloc_align(M * sizeof(reordered_particles), 16, &warray2_mem);
    for(size_t i = 0; i < arr->size; i += M * 4) {
        // load particles into working array:
        for(size_t k = 0; k < M; ++k) {
            size_t base = i + k * 4;
            __m128 wa_p0;
            __m128 wa_p1;
            __m128 wa_p2;
            __m128 wa_m;
            __m128 wa_v0;
            __m128 wa_v1;
            __m128 wa_v2;
            __m128 wa_c;
            for(size_t j = 0; j < 4; ++j) {
                __m128 p_l = _mm_load_ps(&arr->particles[base + j].position[0]);
                __m128 p_h = _mm_load_ps(&arr->particles[base + j].velocity[0]);
                wa_p0 = _mm_insert_ps(wa_p0, p_l, j << 4 | 0 << 6);
                wa_p1 = _mm_insert_ps(wa_p1, p_l, j << 4 | 1 << 6);
                wa_p2 = _mm_insert_ps(wa_p2, p_l, j << 4 | 2 << 6);
                wa_m  = _mm_insert_ps(wa_m , p_l, j << 4 | 3 << 6);
                wa_v0 = _mm_insert_ps(wa_v0, p_h, j << 4 | 0 << 6);
                wa_v1 = _mm_insert_ps(wa_v1, p_h, j << 4 | 1 << 6);
                wa_v2 = _mm_insert_ps(wa_v2, p_h, j << 4 | 2 << 6);
                wa_c  = _mm_insert_ps(wa_c , p_h, j << 4 | 3 << 6);
            }
            _mm_store_ps(&warray1[k].position[0][0], wa_p0);
            _mm_store_ps(&warray1[k].position[1][0], wa_p1);
            _mm_store_ps(&warray1[k].position[2][0], wa_p2);
            _mm_store_ps(&warray1[k].mass    [0][0], wa_m);
            _mm_store_ps(&warray1[k].velocity[0][0], wa_v0);
            _mm_store_ps(&warray1[k].velocity[1][0], wa_v1);
            _mm_store_ps(&warray1[k].velocity[2][0], wa_v2);
            _mm_store_ps(&warray1[k].charge  [0][0], wa_c);
        }
        for(size_t j = 0; effects[j].apply.one != NULL; ++j) {
            if(effects[j].particles == 1) {
                effects[j].apply.one(warray1, effects[j].userdata, dt);
            } else if(effects[j].particles == 2) {
                /*for(size_t k = 0; k<i; ++k) {
                    effects[j].apply.two(&arr->particles[i], &arr->particles[k], effects[j].userdata, dt);
                }*/
            }
        }
        // store particles:
        for(size_t k = 0; k < M; ++k) {
            size_t base = i + k * 4;
            for(size_t j = 0; j < 4; ++j) {
                arr->particles[base + j].position[0] = warray1[k].position[0][j];
                arr->particles[base + j].position[1] = warray1[k].position[1][j];
                arr->particles[base + j].position[2] = warray1[k].position[2][j];
                //arr->particles[base + j].mass        = warray1[k].mass    [0][j]; // cannot change mass
                arr->particles[base + j].velocity[0] = warray1[k].velocity[0][j];
                arr->particles[base + j].velocity[1] = warray1[k].velocity[1][j];
                arr->particles[base + j].velocity[2] = warray1[k].velocity[2][j];
                //arr->particles[base + j].charge      = warray1[k].charge  [0][j]; // cannot change charge
            }
            /*
            size_t base = i + k * 4;
            __m128 wa_p0;
            __m128 wa_p1;
            __m128 wa_p2;
            __m128 wa_m;
            __m128 wa_v0;
            __m128 wa_v1;
            __m128 wa_v2;
            __m128 wa_c;
            wa_p0 = _mm_load_ps(&warray1[k].position[0][0]);
            wa_p1 = _mm_load_ps(&warray1[k].position[1][0]);
            wa_p2 = _mm_load_ps(&warray1[k].position[2][0]);
            wa_m  = _mm_load_ps(&warray1[k].mass    [0][0]);
            wa_v0 = _mm_load_ps(&warray1[k].velocity[0][0]);
            wa_v1 = _mm_load_ps(&warray1[k].velocity[1][0]);
            wa_v2 = _mm_load_ps(&warray1[k].velocity[2][0]);
            wa_c  = _mm_load_ps(&warray1[k].charge  [0][0]);
            for(size_t j = 0; j < 4; ++j) {
                __m128 p_l;
                __m128 p_h;
                p_l = _mm_insert_ps(p_l, wa_p0, j << 6 | 0 << 4);
                p_l = _mm_insert_ps(p_l, wa_p1, j << 6 | 1 << 4);
                p_l = _mm_insert_ps(p_l, wa_p2, j << 6 | 2 << 4);
                p_l = _mm_insert_ps(p_l, wa_m , j << 6 | 3 << 4);
                p_h = _mm_insert_ps(p_h, wa_v0, j << 6 | 0 << 4);
                p_h = _mm_insert_ps(p_h, wa_v1, j << 6 | 1 << 4);
                p_h = _mm_insert_ps(p_h, wa_v2, j << 6 | 2 << 4);
                p_h = _mm_insert_ps(p_h, wa_c , j << 6 | 3 << 4);
                _mm_store_ps(&arr->particles[base + j].position[0], p_l);
                _mm_store_ps(&arr->particles[base + j].velocity[0], p_h);
            }*/
        }
    }
    free(warray1_mem);
    warray1_mem = NULL;
    free(warray2_mem);
    warray2_mem = NULL;
}

void effect_program_c_o3_perf_c(const effect_program *self, const particle_array *arr, float dt, performance_count *out) {
    // set all counters to zero:
    memset(out, 0, sizeof(*out));
    if (self->usr == NULL) {
        // fail, no compilation result.
        return;
    }
    particle_effect_c_o3 *effects = (particle_effect_c_o3 *)self->usr;
    // create working array:
    void *warray1_mem = NULL;
    void *warray2_mem = NULL;
    reordered_particles *warray1 = malloc_align(M * sizeof(reordered_particles), 16, &warray1_mem);
    reordered_particles *warray2 = malloc_align(M * sizeof(reordered_particles), 16, &warray2_mem);
    for(size_t i = 0; i < arr->size; i += M * 4) {
        // load particles into working array:
        for(size_t k = 0; k < M; ++k) {
            size_t base = i + k * 4;
            for(size_t j = 0; j < 4; ++j) {
                warray1[k].position[0][j] = arr->particles[base + j].position[0];
                warray1[k].position[1][j] = arr->particles[base + j].position[1];
                warray1[k].position[2][j] = arr->particles[base + j].position[2];
                warray1[k].mass    [0][j] = arr->particles[base + j].mass;
                warray1[k].velocity[0][j] = arr->particles[base + j].velocity[0];
                warray1[k].velocity[1][j] = arr->particles[base + j].velocity[1];
                warray1[k].velocity[2][j] = arr->particles[base + j].velocity[2];
                warray1[k].charge  [0][j] = arr->particles[base + j].charge;
            }
        }
        for(size_t j = 0; effects[j].apply.one != NULL; ++j) {
            if(effects[j].particles == 1) {
                if (effects[j].perf_c.one != NULL)
                    effects[j].perf_c.one(warray1, effects[j].userdata, dt, out);
            } else if(effects[j].particles == 2) {
                for(size_t k = 0; k<i; ++k) {
                    /*if (effects[j].perf_c.two != NULL)
                        effects[j].perf_c.two(&arr->particles[i], &arr->particles[k], effects[j].userdata, dt, out);*/
                }
            }
        }
    }
    free(warray1_mem);
    warray1_mem = NULL;
    free(warray2_mem);
    warray2_mem = NULL;
}

static void linear_accel_apply(reordered_particles *p, const void *data0, float dt) {
    const float *f_data = (const float *)data0;
    __m128 dtv    = _mm_set_ps1(dt);
    __m128 a0     = _mm_load_ps1(&f_data[0]);
    __m128 a1     = _mm_load_ps1(&f_data[1]);
    __m128 a2     = _mm_load_ps1(&f_data[2]);
    for (size_t k = 0; k < M; ++k) {
        __m128 v0     = _mm_load_ps(&p[k].velocity[0][0]);
        __m128 v1     = _mm_load_ps(&p[k].velocity[1][0]);
        __m128 v2     = _mm_load_ps(&p[k].velocity[2][0]);
        __m128 dt_a0  = _mm_mul_ps(dtv, a0);
        __m128 dt_a1  = _mm_mul_ps(dtv, a1);
        __m128 dt_a2  = _mm_mul_ps(dtv, a2);
        __m128 v1_0   = _mm_add_ps(v0 , dt_a0);
        __m128 v1_1   = _mm_add_ps(v1 , dt_a1);
        __m128 v1_2   = _mm_add_ps(v2 , dt_a2);
        _mm_store_ps(&p[k].velocity[0][0], v1_0);
        _mm_store_ps(&p[k].velocity[1][0], v1_1);
        _mm_store_ps(&p[k].velocity[2][0], v1_2);
    }
}

static void linear_accel_perf_c(const reordered_particles *p, void *data0, float dt, performance_count *out) {
    (void) p; (void) data0; (void) dt;
    out->add += M * 3 * 4;
    out->mul += M * 3 * 4;
    out->loads += 3 + M * 3 * 4;
    out->stores += M * 3 * 4;
}

static particle_effect_c_o3 linear_accel_effect(float x, float y, float z) {
    particle_effect_c_o3 result;
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

static void linear_force_apply(reordered_particles *p, const void *data0, float dt) {
    const float *f_data = (const float *)data0;
    __m128 dtv    = _mm_set_ps1(dt);
    __m128 f0     = _mm_load_ps1(&f_data[0]);
    __m128 f1     = _mm_load_ps1(&f_data[1]);
    __m128 f2     = _mm_load_ps1(&f_data[2]);
    for (size_t k = 0; k < M; ++k) {
        __m128 m      = _mm_load_ps(&p[k].mass    [0][0]);
        __m128 v0     = _mm_load_ps(&p[k].velocity[0][0]);
        __m128 v1     = _mm_load_ps(&p[k].velocity[1][0]);
        __m128 v2     = _mm_load_ps(&p[k].velocity[2][0]);
        __m128 dt_m   = _mm_div_ps(dtv, m);
        __m128 vd0    = _mm_mul_ps(dt_m, f0);
        __m128 vd1    = _mm_mul_ps(dt_m, f1);
        __m128 vd2    = _mm_mul_ps(dt_m, f2);
        __m128 v1_0   = _mm_add_ps(v0 , vd0);
        __m128 v1_1   = _mm_add_ps(v1 , vd1);
        __m128 v1_2   = _mm_add_ps(v2 , vd2);
        _mm_store_ps(&p[k].velocity[0][0], v1_0);
        _mm_store_ps(&p[k].velocity[1][0], v1_1);
        _mm_store_ps(&p[k].velocity[2][0], v1_2);
    }
}

static void linear_force_perf_c(const reordered_particles *p, void *data0, float dt, performance_count *out) {
    (void) p; (void) data0; (void) dt;
    out->add += M * 4 * 3;
    out->mul += M * 4 * 3;
    out->div += M * 4 * 1;
    out->loads += 3 + M * 4 * 4;
    out->stores +=    M * 4 * 3;
}

static particle_effect_c_o3 linear_force_effect(float x, float y, float z) {
    particle_effect_c_o3 result;
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

static void gravitational_force_apply(reordered_particles *p, const void *data0, float dt) {
    const float *f_data = (const float *)data0;
    __m128 dtv    = _mm_set_ps1(dt);
    __m128 cnt0   = _mm_load_ps1(&f_data[0]);
    __m128 cnt1   = _mm_load_ps1(&f_data[1]);
    __m128 cnt2   = _mm_load_ps1(&f_data[2]);
    __m128 mu     = _mm_load_ps1(&f_data[3]);
    __m128 mu_dt  = _mm_mul_ps(mu, dtv);
    for (size_t k = 0; k < M; ++k) {
        __m128 p0     = _mm_load_ps(&p[k].position[0][0]);
        __m128 p1     = _mm_load_ps(&p[k].position[1][0]);
        __m128 p2     = _mm_load_ps(&p[k].position[2][0]);
        __m128 v0     = _mm_load_ps(&p[k].velocity[0][0]);
        __m128 v1     = _mm_load_ps(&p[k].velocity[1][0]);
        __m128 v2     = _mm_load_ps(&p[k].velocity[2][0]);
        __m128 d0     = _mm_sub_ps(p0, cnt0);
        __m128 d1     = _mm_sub_ps(p1, cnt1);
        __m128 d2     = _mm_sub_ps(p2, cnt2);
        __m128 d2_0   = _mm_mul_ps(d0, d0);
        __m128 d2_1   = _mm_mul_ps(d1, d1);
        __m128 d2_2   = _mm_mul_ps(d2, d2);
        __m128 r2_1   = _mm_add_ps(d2_0, d2_1);
        __m128 r2     = _mm_add_ps(r2_1, d2_2);
        __m128 d0_mu_dt=_mm_mul_ps(d0, mu_dt);
        __m128 d1_mu_dt=_mm_mul_ps(d1, mu_dt);
        __m128 d2_mu_dt=_mm_mul_ps(d2, mu_dt);
        __m128 r      = _mm_sqrt_ps(r2);
        __m128 r4     = _mm_mul_ps(r2, r2);
        // reciprocal with one newton step for accuracy:
        __m128 r4_ra  = _mm_rcp_ps(r4);
        __m128 r4_r2  = _mm_add_ps(r4_ra, r4_ra);
        __m128 r8_ra  = _mm_mul_ps(r4_ra, r4_ra);
        __m128 r4_r0  = _mm_mul_ps(r8_ra, r4);
        __m128 r4_r   = _mm_sub_ps(r4_r2, r4_r0);
        
        __m128 d0_r   = _mm_mul_ps(d0_mu_dt, r);
        __m128 d1_r   = _mm_mul_ps(d1_mu_dt, r);
        __m128 d2_r   = _mm_mul_ps(d2_mu_dt, r);
        
        __m128 d0_k_2 = _mm_mul_ps(d0_r, r4_r);
        __m128 d1_k_2 = _mm_mul_ps(d1_r, r4_r);
        __m128 d2_k_2 = _mm_mul_ps(d2_r, r4_r);
        
        __m128 v1_0   = _mm_add_ps(v0, d0_k_2);
        __m128 v1_1   = _mm_add_ps(v1, d1_k_2);
        __m128 v1_2   = _mm_add_ps(v2, d2_k_2);
        
        _mm_store_ps(&p[k].velocity[0][0], v1_0);
        _mm_store_ps(&p[k].velocity[1][0], v1_1);
        _mm_store_ps(&p[k].velocity[2][0], v1_2);
    }
}

static void gravitational_force_perf_c(const reordered_particles *p, void *data0, float dt, performance_count *out) {
    (void) p; (void) data0; (void) dt;
    out->add +=    M * 4 * 10;
    out->mul +=    M * 4 * 15 + 1;
    out->rcp +=    M * 4 * 1;
    out->sqrt +=   M * 4 * 1;
    out->loads +=  M * 4 * 6 + 4;
    out->stores += M * 4 * 3;
}

static particle_effect_c_o3 gravitational_force_effect(float x, float y, float z, float mu) {
    particle_effect_c_o3 result;
    result.particles = 1;
    result.apply.one  = gravitational_force_apply;
    result.perf_c.one = gravitational_force_perf_c;
    float *data = malloc_align(4*sizeof(float), 16, &result.p_to_free);
    data[0] = x;
    data[1] = y;
    data[2] = z;
    data[3] = mu;
    result.userdata = data;
    result.p_to_free = data;
    return result;
}

static void plane_bounce_apply(reordered_particles *p, const void *data0, float dt) {
    (void) dt;
    const float *f_data = (const float *)data0;
    __m128 nv0    = _mm_load_ps1(&f_data[0]);
    __m128 nv1    = _mm_load_ps1(&f_data[1]);
    __m128 nv2    = _mm_load_ps1(&f_data[2]);
    __m128 d      = _mm_load_ps1(&f_data[3]);
    __m128 dec    = _mm_load_ps1(&f_data[4]);
    
    __m128 nv0_2  = _mm_mul_ps(nv0, _mm_set_ps1(2.0f));
    __m128 nv1_2  = _mm_mul_ps(nv1, _mm_set_ps1(2.0f));
    __m128 nv2_2  = _mm_mul_ps(nv2, _mm_set_ps1(2.0f));
    
    __m128 z      = _mm_setzero_ps();
    __m128 one    = _mm_set_ps1(1.0f);
    for (size_t k = 0; k < M; ++k) {
        __m128 p0     = _mm_load_ps(&p[k].position[0][0]);
        __m128 p1     = _mm_load_ps(&p[k].position[1][0]);
        __m128 p2     = _mm_load_ps(&p[k].position[2][0]);
        __m128 v0     = _mm_load_ps(&p[k].velocity[0][0]);
        __m128 v1     = _mm_load_ps(&p[k].velocity[1][0]);
        __m128 v2     = _mm_load_ps(&p[k].velocity[2][0]);
        
        __m128 pd0    = _mm_mul_ps(p0, nv0);
        __m128 pd1    = _mm_mul_ps(p1, nv1);
        __m128 pd2    = _mm_mul_ps(p2, nv2);
        __m128 vd0    = _mm_mul_ps(v0, nv0);
        __m128 vd1    = _mm_mul_ps(v1, nv1);
        __m128 vd2    = _mm_mul_ps(v2, nv2);
        
        __m128 ddst_0 = _mm_add_ps(pd0, pd1);
        __m128 ddst_1 = _mm_sub_ps(pd2, d);
        __m128 vnrm_0 = _mm_add_ps(vd0, vd1);
        __m128 ddst   = _mm_add_ps(ddst_0, ddst_1);
        __m128 vnrm   = _mm_add_ps(vnrm_0, vd2);
        
        __m128 ddst_z = _mm_min_ps(z, ddst);
        __m128 vnrm_z = _mm_min_ps(z, vnrm);
        
        __m128 vmask  = _mm_cmplt_ps(ddst, z);  // mask is true when particle should be modified
        __m128 dmask  = _mm_cmplt_ps(vnrm, z);  // mask is true when particle should be modified
        
        __m128 ddst_zm= _mm_and_ps(ddst_z, dmask);
        __m128 vnrm_zm= _mm_and_ps(vnrm_z, vmask);
        
        __m128 ddst_n0= _mm_mul_ps(ddst_zm, nv0);
        __m128 ddst_n1= _mm_mul_ps(ddst_zm, nv1);
        __m128 ddst_n2= _mm_mul_ps(ddst_zm, nv2);
        
        __m128 dec_m1 = _mm_blendv_ps(one, dec,    vmask);
        __m128 dec_m  = _mm_blendv_ps(one, dec_m1, dmask);
        
        __m128 vnrm_2n0= _mm_mul_ps(vnrm_zm, nv0_2);
        __m128 vnrm_2n1= _mm_mul_ps(vnrm_zm, nv1_2);
        __m128 vnrm_2n2= _mm_mul_ps(vnrm_zm, nv2_2);
        
        __m128 v1_0   = _mm_sub_ps(v0, vnrm_2n0);
        __m128 v1_1   = _mm_sub_ps(v1, vnrm_2n1);
        __m128 v1_2   = _mm_sub_ps(v2, vnrm_2n2);
        
        __m128 p1_0   = _mm_sub_ps(p0, ddst_n0);
        __m128 p1_1   = _mm_sub_ps(p1, ddst_n1);
        __m128 p1_2   = _mm_sub_ps(p2, ddst_n2);
        
        __m128 v2_0   = _mm_mul_ps(v1_0, dec_m);
        __m128 v2_1   = _mm_mul_ps(v1_1, dec_m);
        __m128 v2_2   = _mm_mul_ps(v1_2, dec_m);
        
        _mm_store_ps(&p[k].position[0][0], p1_0);
        _mm_store_ps(&p[k].position[1][0], p1_1);
        _mm_store_ps(&p[k].position[2][0], p1_2);
        _mm_store_ps(&p[k].velocity[0][0], v2_0);
        _mm_store_ps(&p[k].velocity[1][0], v2_1);
        _mm_store_ps(&p[k].velocity[2][0], v2_2);
    }
}

static void plane_bounce_perf_c(const particle *p, void *data0, float dt, performance_count *out) {
    (void) dt;
    float *data = data0;
    float dist =    data[0]*p->position[0] + data[1]*p->position[1] + data[2]*p->position[2];
    float vnormal = data[0]*p->velocity[0] + data[1]*p->velocity[1] + data[2]*p->velocity[2];
    float d = data[3];
    
    out->add +=    M * 4 * 11;
    out->mul +=    M * 4 * 15 + 3 * 4;
    out->cmp +=    M * 4 * 4;
    out->loads +=  M * 4 * 6 + 5;
    out->stores += M * 4 * 6;
}

static particle_effect_c_o3 plane_bounce_effect(float x, float y, float z, float d, float a) {
    particle_effect_c_o3 result;
    result.particles = 1;
    result.apply.one =  plane_bounce_apply;
    result.perf_c.one = plane_bounce_perf_c;
    float *data = malloc_align(5*sizeof(float), 16, &result.p_to_free);
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

    //float fd0 = f_data[0];
    //float fd1 = f_data[1];
    //float fd2 = f_data[2];
    //float pp0 = p->position[0];
    //float pp1 = p->position[1];
    //float pp2 = p->position[2];
    //float pv0 = p->velocity[0];
    //float pv1 = p->velocity[1];
    //float pv2 = p->velocity[2];
    __m128 z        = _mm_set_ps1(0.);
    __m128 ppt      = _mm_load_ps(&p->position[0]);
    __m128 pp       = _mm_blend_ps(ppt, z, 0b1000);
    __m128 pvt      = _mm_load_ps(&p->velocity[0]);
    __m128 pv       = _mm_blend_ps(pvt, z, 0b1000);
    __m128 fdt      = _mm_load_ps(&f_data[0]);
    __m128 fd       = _mm_blend_ps(fdt, z, 0b1000);

    float d   = f_data[3];
    float dec = f_data[4];
    float d2 = d * d;
    //float n0  = pp0-fd0;
    //float n1  = pp1-fd1;
    //float n2  = pp2-fd2;
    __m128 n        = _mm_sub_ps(pp, fd);

    //float r2 = n0*n0 + n1*n1 + n2*n2;
    __m128 r2m      = _mm_mul_ps(n,n);
    __m128 r2t      = _mm_hadd_ps(r2m,r2m);
    __m128 r2       = _mm_hadd_ps(r2t,r2t);

    //float vnormal = n0*pv0 + n1*pv1 + n2*pv2;
    __m128 vntt     = _mm_mul_ps(pv, n);
    __m128 vnt      = _mm_hadd_ps(vntt,vntt);
    __m128 vn       = _mm_hadd_ps(vnt,vnt);

    //float r2_reci = 1.0f / r2;
    __m128 r2_reci  = _mm_div_ps(_mm_set_ps1(1.0),r2);
    float r2s, vns;
    _mm_store_ss(&r2s,r2);
    _mm_store_ss(&vns, vn);
    if(d2 > r2s && vns < 0.0f) {
        // inside sphere and going further inside sphere:
        //float r_reci = sqrtf(r2_reci);
        __m128 r_reci   = _mm_sqrt_ps(r2_reci);

        //float vnormal2_r2 = 2.0f * vnormal * r2_reci;
        __m128 vn2_r2t  = _mm_mul_ps(_mm_set_ps1(2.0), vn);
        __m128 vn2_r2   = _mm_mul_ps(vn2_r2t, r2_reci);

        //float pv0_1 = pv0 - vnormal2_r2*n0;
        //float pv1_1 = pv1 - vnormal2_r2*n1;
        //float pv2_1 = pv2 - vnormal2_r2*n2;

        __m128 vn_n     = _mm_mul_ps(vn2_r2, n);
        __m128 pvut     = _mm_sub_ps(pv, vn_n);

        //float d_r_r_reci = d * r_reci - 1.0f;
        __m128 dt           = _mm_load_ps1(&d);
        __m128 d_reci_1t    = _mm_mul_ps(dt,r_reci);
        __m128 d_reci_1     = _mm_sub_ps(d_reci_1t,_mm_set_ps1(1.0));

        //p->velocity[0] = pv0_1 * dec;
        //p->velocity[1] = pv1_1 * dec;
        //p->velocity[2] = pv2_1 * dec;
        __m128 vdec     = _mm_load_ps1(&dec);
        __m128 pvu      = _mm_mul_ps(vdec, pvut);
        __m128 pv4      = _mm_load_ps1(&p->velocity[3]);
        __m128 pvus     = _mm_blend_ps(pvu, pv4, 0b1000);
        _mm_store_ps(&p->velocity[0], pvus);

        //p->position[0] = pp0 + d_r_r_reci*n0;
        //p->position[1] = pp1 + d_r_r_reci*n1;
        //p->position[2] = pp2 + d_r_r_reci*n2;
        __m128 d_nt     = _mm_mul_ps(n, d_reci_1);
        __m128 d_n      = _mm_add_ps(pp, d_nt);
        __m128 pp4      = _mm_load_ps1(&p->position[3]);
        __m128 ppu      = _mm_blend_ps(d_n, pp4, 0b1000);
        _mm_store_ps(&p->position[0], ppu);
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

static particle_effect_c_o3 sphere_bounce_effect(float x, float y, float z, float r, float a) {
    particle_effect_c_o3 result;
    result.particles = 1;
    result.apply.one =  sphere_bounce_apply;
    result.perf_c.one = sphere_bounce_perf_c;
    float *data = malloc_align(5*sizeof(float), 16, &result.p_to_free);
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

    __m128 z        = _mm_set_ps1(0.);
    __m128 pp1      = _mm_load_ps(&p1->position[0]);
    __m128 pp2      = _mm_load_ps(&p2->position[0]);
    __m128 pv1      = _mm_load_ps(&p1->velocity[0]);
    __m128 pv2      = _mm_load_ps(&p2->velocity[0]);

    float dtmu = dt * mu;
    __m128 m1       = _mm_load_ps1(&p1->mass);
    __m128 m2       = _mm_load_ps1(&p2->mass);
    __m128 dt_mu    = _mm_load_ps1(&dtmu);

    //float d0 = pp20-pp10;
    //float d1 = pp21-pp11;
    //float d2 = pp22-pp12;
    //float r2 = d0*d0 + d1*d1 + d2*d2;
    //float r4 = r2 * r2;
    //float r = sqrtf(r2);
    __m128 dift     = _mm_sub_ps(pp2,pp1);
    __m128 dif      = _mm_blend_ps(dift, z, 0b1000);
    __m128 dif2     = _mm_mul_ps(dif,dif);
    __m128 r2t      = _mm_hadd_ps(dif2,dif2);
    __m128 r2       = _mm_hadd_ps(r2t,r2t);
    __m128 r4       = _mm_mul_ps(r2,r2);
    __m128 r        = _mm_sqrt_ps(r2);
    // is r^3 because we need to normalize diff.
    //float mu_0 = dt_mu / r4;
    //float mu_1 = mu_0 * r;
    //float mu_m2 = mu_1 * m2;
    //float mu_m1 = mu_1 * m1;

    __m128 mu_0     = _mm_div_ps(dt_mu, r4);
    __m128 mu_1     = _mm_mul_ps(mu_0, r);
    __m128 mu_m1    = _mm_mul_ps(mu_1,m1);
    __m128 mu_m2    = _mm_mul_ps(mu_1,m2);

    //float d0_1 = d0 * mu_m2;
    //float d1_1 = d1 * mu_m2;
    //float d2_1 = d2 * mu_m2;
    //float d0_2 = d0 * mu_m1;
    //float d1_2 = d1 * mu_m1;
    //float d2_2 = d2 * mu_m1;

    __m128 dv1t     = _mm_mul_ps(dif, mu_m2);
    __m128 dv2t     = _mm_mul_ps(dif, mu_m1);

    //p1->velocity[0] = pv10 - d0_1;
    //p1->velocity[1] = pv11 - d1_1;
    //p1->velocity[2] = pv12 - d2_1;
    //p2->velocity[0] = pv20 + d0_2;
    //p2->velocity[1] = pv21 + d1_2;
    //p2->velocity[2] = pv22 + d2_2;
    __m128 dv1u     = _mm_sub_ps(pv1,dv1t);
    __m128 dv2u     = _mm_add_ps(pv2,dv2t);

    __m128 c1       = _mm_load_ps1(&p1->charge);
    __m128 c2       = _mm_load_ps1(&p2->charge);
    __m128 dv1      = _mm_blend_ps(dv1u, c1, 0b1000);
    __m128 dv2      = _mm_blend_ps(dv2u, c2, 0b1000);

    _mm_store_ps(&p1->velocity[0], dv1);
    _mm_store_ps(&p2->velocity[0], dv2);


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

static particle_effect_c_o3 pairwise_gravitational_force_effect(float mu) {
    particle_effect_c_o3 result;
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
    __m128 z        = _mm_set_ps1(0.);
    __m128 pp1      = _mm_load_ps(&p1->position[0]);
    __m128 pp2      = _mm_load_ps(&p2->position[0]);
    __m128 pv1      = _mm_load_ps(&p1->velocity[0]);
    __m128 pv2      = _mm_load_ps(&p2->velocity[0]);
    //float pp10 = p1->position[0];
    //float pp11 = p1->position[1];
    //float pp12 = p1->position[2];
    //float pp20 = p2->position[0];
    //float pp21 = p2->position[1];
    //float pp22 = p2->position[2];
    //float pv10 = p1->velocity[0];
    //float pv11 = p1->velocity[1];
    //float pv12 = p1->velocity[2];
    //float pv20 = p2->velocity[0];
    //float pv21 = p2->velocity[1];
    //float pv22 = p2->velocity[2];
    //float d0 = pp20 - pp10;
    //float d1 = pp21 - pp11;
    //float d2 = pp22 - pp12;
    //float s_radius2 = radius * radius;
    //float s_r2 = d0 * d0 + d1 * d1 + d2 * d2;

    __m128 dift     = _mm_sub_ps(pp2, pp1);
    __m128 difp     = _mm_blend_ps(dift, z, 0b1000); //[d0 d1 d2 0]
    __m128 difp2    = _mm_mul_ps(difp,difp);
    __m128 difp2ss  = _mm_hadd_ps(difp2,difp2);
    __m128 sump2s   = _mm_hadd_ps(difp2ss, difp2ss); // sum^2

    // project velocities onto diff vector:
    //float pv10d0 = pv10*d0;
    //float pv11d1 = pv11*d1;
    //float pv12d2 = pv12*d2;
    //float pv20d0 = pv20*d0;
    //float pv21d1 = pv21*d1;
    //float pv22d2 = pv22*d2;
    __m128 propv1   = _mm_mul_ps(pv1, difp);
    __m128 propv2   = _mm_mul_ps(pv2, difp);

    //float u1_0 = pv10d0 + pv11d1 + pv12d2;
    //float u2_0 = pv20d0 + pv21d1 + pv22d2;
    __m128 u_0t     = _mm_hadd_ps(propv1,propv2);
    __m128 u_0      = _mm_hadd_ps(u_0t, u_0t); // u0 = [u1_0, u2_0, u1_0, u2_0]

    float u0[4];
    _mm_store_ps(&u0[0], u_0);

    float radius2 = radius * radius;
    float r2;
    _mm_store_ss(&r2, sump2s);
    //float u1_u2_0 = u1_0 - u2_0;
    float u1_u2_0 = u0[0] - u0[1];

    if(r2 < radius2) {
        // within collision range:
        float r = sqrtf(r2);
        float r_reci = 1.0f / r;
        float radius_r = radius - r;
        float radius_r05 = 0.5f * radius_r;

        //float d0_1 = d0 * r_reci;
        //float d1_1 = d1 * r_reci;
        //float d2_1 = d2 * r_reci;
        __m128 r_rec    = _mm_load_ps1(&r_reci);
        __m128 dif1     = _mm_mul_ps(difp, r_rec);

        if(u1_u2_0 > 0.0f) {
            float m1 = p1->mass;
            float m2 = p2->mass;
            float restitution = f_data[1];
            float m12_reci = 1.0f/(m1+m2);
            float u1 = u0[0] * r_reci;
            float u2 = u0[1] * r_reci;
            float um1pum2 = u1*m1+u2*m2;
            float u1_u2 = u1_u2_0 * r_reci;
            // particles are approaching:
            // calculate new velocities according to collision laws: (in the direction of diff)
            float v1 = (um1pum2-restitution*m2*u1_u2)*m12_reci;
            float v2 = (um1pum2+restitution*m1*u1_u2)*m12_reci;
            // update velocities:
            float v1_u1 = v1 - u1;
            float v2_u2 = v2 - u2;

            __m128 v1u1     = _mm_load_ps1(&v1_u1);
            __m128 v2u2     = _mm_load_ps1(&v2_u2);

            //float d0_1_v1_u1 = d0_1 * v1_u1;
            //float d1_1_v1_u1 = d1_1 * v1_u1;
            //float d2_1_v1_u1 = d2_1 * v1_u1;
            //float d0_1_v2_u2 = d0_1 * v2_u2;
            //float d1_1_v2_u2 = d1_1 * v2_u2;
            //float d2_1_v2_u2 = d2_1 * v2_u2;

            __m128 d1v1u1   = _mm_mul_ps(dif1, v1u1);
            __m128 d1v2u2   = _mm_mul_ps(dif1, v2u2);

            //p1->velocity[0] = pv10 + d0_1_v1_u1;
            //p1->velocity[1] = pv11 + d1_1_v1_u1;
            //p1->velocity[2] = pv12 + d2_1_v1_u1;
            //p2->velocity[0] = pv20 + d0_1_v2_u2;
            //p2->velocity[1] = pv21 + d1_1_v2_u2;
            //p2->velocity[2] = pv22 + d2_1_v2_u2;
            __m128 v1uv     = _mm_add_ps(pv1, d1v1u1);
            __m128 v2uv     = _mm_add_ps(pv2, d1v2u2);
            __m128 v1_rec   = _mm_load_ps1(&p1->velocity[3]);
            __m128 v2_rec   = _mm_load_ps1(&p2->velocity[3]);
            __m128 v1_rec1  = _mm_blend_ps(v1uv, v1_rec, 0b1000);
            __m128 v2_rec1  = _mm_blend_ps(v2uv, v2_rec, 0b1000);

            _mm_store_ps(&p1->velocity[0], v1_rec1);
            _mm_store_ps(&p2->velocity[0], v2_rec1);
        }
        // move particles apart, such that they don't intersect (in the direction of diff)
        //float d0_1_rr05 = d0_1*radius_r05;
        //float d1_1_rr05 = d1_1*radius_r05;
        //float d2_1_rr05 = d2_1*radius_r05;

        //p1->position[0] = pp10 - d0_1_rr05;
        //p1->position[1] = pp11 - d1_1_rr05;
        //p1->position[2] = pp12 - d2_1_rr05;
        //p2->position[0] = pp20 + d0_1_rr05;
        //p2->position[1] = pp21 + d1_1_rr05;
        //p2->position[2] = pp22 + d2_1_rr05;

        __m128 r05      = _mm_load_ps1(&radius_r05);
        __m128 dr05     = _mm_mul_ps(dif1, r05);

        __m128 pp1_u    = _mm_sub_ps(pp1, dr05);
        __m128 pp2_u    = _mm_add_ps(pp2, dr05);


        __m128 p1_rec   = _mm_load_ps1(&p1->position[3]);
        __m128 p2_rec   = _mm_load_ps1(&p2->position[3]);
        __m128 p1_rec1  = _mm_blend_ps(pp1_u, p1_rec, 0b1000);
        __m128 p2_rec1  = _mm_blend_ps(pp2_u, p2_rec, 0b1000);
        _mm_store_ps(&p1->position[0], p1_rec1);
        _mm_store_ps(&p2->position[0], p2_rec1);
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

static particle_effect_c_o3 pairwise_sphere_collision_effect(float radius, float restitution) {
    particle_effect_c_o3 result;
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

static void newton_step_apply(reordered_particles *p, const void *data0, float dt) {
    (void) data0;
    const float *f_data = (const float *)data0;
    __m128 dtv    = _mm_set_ps1(dt);
    for (size_t k = 0; k < M; ++k) {
        __m128 v0     = _mm_load_ps(&p[k].velocity[0][0]);
        __m128 v1     = _mm_load_ps(&p[k].velocity[1][0]);
        __m128 v2     = _mm_load_ps(&p[k].velocity[2][0]);
        __m128 p0     = _mm_load_ps(&p[k].position[0][0]);
        __m128 p1     = _mm_load_ps(&p[k].position[1][0]);
        __m128 p2     = _mm_load_ps(&p[k].position[2][0]);
        __m128 dt_v0  = _mm_mul_ps(dtv, v0);
        __m128 dt_v1  = _mm_mul_ps(dtv, v1);
        __m128 dt_v2  = _mm_mul_ps(dtv, v2);
        __m128 p1_0   = _mm_add_ps(p0 , dt_v0);
        __m128 p1_1   = _mm_add_ps(p1 , dt_v1);
        __m128 p1_2   = _mm_add_ps(p2 , dt_v2);
        _mm_store_ps(&p[k].position[0][0], p1_0);
        _mm_store_ps(&p[k].position[1][0], p1_1);
        _mm_store_ps(&p[k].position[2][0], p1_2);
    }
}

static void newton_step_perf_c(const reordered_particles *p, void *data0, float dt, performance_count *out) {
    (void) p; (void) data0; (void) dt;
    out->add +=    M * 4 * 3;
    out->mul +=    M * 4 * 3;
    out->loads +=  M * 4 * 6;
    out->stores += M * 4 * 3;
}

static particle_effect_c_o3 newton_step_effect() {
    particle_effect_c_o3 result;
    result.particles = 1;
    result.apply.one  = newton_step_apply;
    result.perf_c.one = newton_step_perf_c;
    result.userdata = NULL;
    result.p_to_free = NULL;
    return result;
}

