#include <particle_effects.h>
#include <math.h>

void linear_force_apply(particle *p, void *data0, float dt) {
    float *data = data0;
    p->velocity[0] += dt*data[0];
    p->velocity[1] += dt*data[1];
    p->velocity[2] += dt*data[2];
}

particle_effect linear_force_effect(float x, float y, float z) {
    particle_effect result;
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

particle_effect gravitational_force_effect(float x, float y, float z, float mu) {
    particle_effect result;
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
    float *data = data0;
    float dist = data[0]*p->position[0] + data[1]*p->position[1] + data[2]*p->position[2];
    float vnormal = data[0]*p->velocity[0] + data[1]*p->velocity[1] + data[2]*p->velocity[2];
    float d = data[3];
    if(dist<d && vnormal<0) {
        p->velocity[0] -= 2.0*vnormal*data[0];
        p->velocity[1] -= 2.0*vnormal*data[1];
        p->velocity[2] -= 2.0*vnormal*data[2];
    }
}

particle_effect plane_bounce_effect(float x, float y, float z, float d) {
    particle_effect result;
    result.apply = plane_bounce_apply;
    float *data = malloc(4*sizeof(float));
    float r = sqrt(x*x + y*y + z*z);
    data[0] = x/r;
    data[1] = y/r;
    data[2] = z/r;
    data[3] = d/r;
    result.userdata = data;
    return result;
}

void sphere_bounce_apply(particle *p, void *data0, float dt) {
    float *data = data0;
    float normal[3];
    normal[0] = p->position[0]-data[0];
    normal[1] = p->position[1]-data[1];
    normal[2] = p->position[2]-data[2];
    double r = normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2];
    normal[0] /= r;
    normal[1] /= r;
    normal[2] /= r;
    float vnormal = normal[0]*p->velocity[0] + normal[1]*p->velocity[1] + normal[2]*p->velocity[2];
    float d = data[3];
    if(r<d && vnormal<0) {
        p->velocity[0] -= 2.0*vnormal*normal[0];
        p->velocity[1] -= 2.0*vnormal*normal[1];
        p->velocity[2] -= 2.0*vnormal*normal[2];
    }
}

particle_effect sphere_bounce_effect(float x, float y, float z, float r) {
    particle_effect result;
    result.apply = sphere_bounce_apply;
    float *data = malloc(4*sizeof(float));
    data[0] = x;
    data[1] = y;
    data[2] = z;
    data[3] = r;
    result.userdata = data;
    return result;
}

void newton_step_apply(particle *p, void *data0, float dt) {
    (void) data0;
    p->position[0] += dt*p->velocity[0];
    p->position[1] += dt*p->velocity[1];
    p->position[2] += dt*p->velocity[2];
}

particle_effect newton_step_effect() {
    particle_effect result;
    result.apply = newton_step_apply;
    result.userdata = NULL;
    return result;

}

