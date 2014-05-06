#ifndef __effect_desc_H
#define __effect_desc_H

#include <stdint.h>
#include <stddef.h>

typedef enum effect_type_t {
    EFFECT_TYPE_LINEAR_ACCEL,   ///< Linear acceleration effect
    EFFECT_TYPE_LINEAR_FORCE,   ///< Linear force effect
    EFFECT_TYPE_CENTRAL_FORCE,  ///< Central force effect
    EFFECT_TYPE_PLANE_BOUNCE,   ///< Plane bounce effect
    EFFECT_TYPE_SPHERE_BOUNCE,  ///< Sphere bounce effect
    EFFECT_TYPE_NEWTON_STEP,    ///< Newton simulation step effect
    EFFECT_TYPE_GRAVITY_FORCE,  ///< pairwise gravitational force
    EFFECT_TYPE_SPHERE_COLLISION,///< particle collision

    EFFECT_TYPE_COUNT
} effect_type;

/** Effect description element struct
 */
typedef struct effect_desc_ele_t {
    effect_type type;       ///< Effect type
    float float_usr[8];     ///< Float user data of the effect.
    void *usr;              ///< Arbitrary user data of the effect.
} effect_desc_ele;

/** Effect description array context
 */
typedef struct effect_desc_t {
    size_t size;                ///< Number of elements in the element array
    size_t capacity;            ///< Capacity of the element array
    effect_desc_ele *elements;  ///< The element array
} effect_desc;

/** Initializes the effect description array context
 *  @return 0 on success, non-zero on error.
 */
int effect_desc_init(effect_desc *ctx);
/** Destroys the effect description array context
 */
void effect_desc_destroy(effect_desc *ctx);
/** Reserves space in the effect description array context
 *  @return 0 on success, non-zero on error.
 */
int effect_desc_reserve(effect_desc *ctx, size_t capacity);
size_t effect_desc_size(const effect_desc *ctx);

/** Removes an effect descriptor element.
 *  @param ctx Effect description array context
 *  @param idx Element index to remove.
 *  @return 0 on success, non-zero on error.
 */
int effect_desc_remove(effect_desc *ctx, int idx);

/** Adds an linear acceleration effect.
 *  @param ctx Effect description array context
 *  @param x, y, z Linear acceleration vector.
 *  @return effect index on success, -1 on error.
 */
int effect_desc_add_linear_accel (effect_desc *ctx, float x, float y, float z);
/** Adds an linear force effect.
 *  @param ctx Effect description array context
 *  @param x, y, z Linear force vector.
 *  @return effect index on success, -1 on error.
 */
int effect_desc_add_linear_force (effect_desc *ctx, float x, float y, float z);
/** Adds an central force effect.
 *  The force increase with the square of the distance.
 *  @param ctx Effect description array context
 *  @param x, y, z Central point of the force.
 *  @param mu Force factor.
 *  @return effect index on success, -1 on error.
 */
int effect_desc_add_central_force(effect_desc *ctx, float x, float y, float z, float mu);
/** Adds a plane bounce effect.
 *  The particles will bounce off the the plane.
 *  @param ctx Effect description array context
 *  @param x, y, z Plane normal vector.
 *  @param d Plane offset.
 *  @param a Bounce attenuation. 1.0 = no attenuation.
 *  @return effect index on success, -1 on error.
 */
int effect_desc_add_plane_bounce (effect_desc *ctx, float x, float y, float z, float d, float a);
/** Adds a sphere bounce effect.
 *  The particles will bounce off the the sphere.
 *  @param ctx Effect description array context
 *  @param x, y, z Sphere centre point
 *  @param r Sphere diameter
 *  @param a Bounce attenuation. 1.0 = no attenuation.
 *  @return effect index on success, -1 on error.
 */
int effect_desc_add_sphere_bounce(effect_desc *ctx, float x, float y, float z, float r, float a);
/** Adds a pairwise gravitation force between particles.
 *  Applies a acceleration of mu*m1*m1/r^2
 *  @param ctx Effect description array context
 *  @param mu force constant
 *  @return effect index on success, -1 on error.
 */
int effect_desc_add_gravity_force(effect_desc *ctx, float mu);
/** Adds a newton simulation step effect.
 *  This will update the particle positions according to their velocity.
 *  @param ctx Effect description array context
 *  @return effect index on success, -1 on error.
 */
int effect_desc_add_sphere_collision(effect_desc *ctx, float radius, float restitution);
/** Adds a newton simulation step effect.
 *  This will update the particle positions according to their velocity.
 *  @param ctx Effect description array context
 *  @return effect index on success, -1 on error.
 */
int effect_desc_add_newton_step  (effect_desc *ctx);




#endif // __effect_desc_H
