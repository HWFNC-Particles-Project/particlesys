#include "effect_desc.h"
#include <stdlib.h>
#include <string.h>

int effect_desc_add_element(effect_desc *ctx, effect_desc_ele *e);

int effect_desc_init(effect_desc *ctx) {
    ctx->size = 0;
    ctx->capacity = 0;
    ctx->elements = NULL;
    return effect_desc_reserve(ctx, 64);
}

void effect_desc_destroy(effect_desc *ctx) {
    // remove all elements:
    while(ctx->size > 0) {
        effect_desc_remove(ctx, ctx->size - 1);
    }
    // free memory:
    free(ctx->elements);
    ctx->elements = NULL;
}

size_t effect_desc_size(const effect_desc *ctx) {
    return ctx->size;
}

int effect_desc_reserve(effect_desc *ctx, size_t capacity) {
    if(capacity < ctx->capacity) {
        return 0;
    }
    effect_desc_ele *new_ele = realloc(ctx->elements, capacity * sizeof(effect_desc_ele));
    if(new_ele != NULL) {
        ctx->capacity = capacity;
        ctx->elements = new_ele;
        return 0;
    } else {
        return 1;
    }
}

int effect_desc_add_element(effect_desc *ctx, effect_desc_ele *e) {
    if(ctx->size >= ctx->capacity) {
        if(effect_desc_reserve(ctx, ctx->capacity/2*3)) {
            return -1;
        }
    }
    int idx = ctx->size;
    ctx->size++;
    ctx->elements[idx] = *e;
    return idx;
}

int effect_desc_remove(effect_desc *ctx, int idx) {
    if (idx < 0 || idx >= ctx->size) return -1;
    effect_desc_ele *e = &ctx->elements[idx];
    switch(e->type) {
        case EFFECT_TYPE_LINEAR_ACCEL:
        case EFFECT_TYPE_LINEAR_FORCE:
        case EFFECT_TYPE_CENTRAL_FORCE:
        case EFFECT_TYPE_PLANE_BOUNCE:
        case EFFECT_TYPE_SPHERE_BOUNCE:
        case EFFECT_TYPE_NEWTON_STEP:
            // nothing todo for this effect type.
            break;
    }
    ctx->size--;
    return 0;
}

int effect_desc_add_linear_accel (effect_desc *ctx, float x, float y, float z) {
    effect_desc_ele e;
    memset(&e, 0, sizeof(effect_desc_ele));
    e.type = EFFECT_TYPE_LINEAR_ACCEL;
    e.float_usr[0] = x;
    e.float_usr[1] = y;
    e.float_usr[2] = z;
    return effect_desc_add_element(ctx, &e);
}

int effect_desc_add_linear_force (effect_desc *ctx, float x, float y, float z) {
    effect_desc_ele e;
    memset(&e, 0, sizeof(effect_desc_ele));
    e.type = EFFECT_TYPE_LINEAR_FORCE;
    e.float_usr[0] = x;
    e.float_usr[1] = y;
    e.float_usr[2] = z;
    return effect_desc_add_element(ctx, &e);
}

int effect_desc_add_central_force(effect_desc *ctx, float x, float y, float z, float mu) {
    effect_desc_ele e;
    memset(&e, 0, sizeof(effect_desc_ele));
    e.type = EFFECT_TYPE_CENTRAL_FORCE;
    e.float_usr[0] = x;
    e.float_usr[1] = y;
    e.float_usr[2] = z;
    e.float_usr[3] = mu;
    return effect_desc_add_element(ctx, &e);
}

int effect_desc_add_plane_bounce (effect_desc *ctx, float x, float y, float z, float d, float a) {
    effect_desc_ele e;
    memset(&e, 0, sizeof(effect_desc_ele));
    e.type = EFFECT_TYPE_PLANE_BOUNCE;
    e.float_usr[0] = x;
    e.float_usr[1] = y;
    e.float_usr[2] = z;
    e.float_usr[3] = d;
    e.float_usr[4] = a;
    return effect_desc_add_element(ctx, &e);
}

int effect_desc_add_sphere_bounce(effect_desc *ctx, float x, float y, float z, float r, float a) {
    effect_desc_ele e;
    memset(&e, 0, sizeof(effect_desc_ele));
    e.type = EFFECT_TYPE_SPHERE_BOUNCE;
    e.float_usr[0] = x;
    e.float_usr[1] = y;
    e.float_usr[2] = z;
    e.float_usr[3] = r;
    e.float_usr[4] = a;
    return effect_desc_add_element(ctx, &e);
}

int effect_desc_add_newton_step  (effect_desc *ctx) {
    effect_desc_ele e;
    memset(&e, 0, sizeof(effect_desc_ele));
    e.type = EFFECT_TYPE_NEWTON_STEP;
    return effect_desc_add_element(ctx, &e);
}

