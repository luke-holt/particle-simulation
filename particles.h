#ifndef PARTICLES_H
#define PARTICLES_H

#include <stddef.h>

#define DIMS (2)

typedef struct {
    float gravity; // system gravity (m/s^2)
    float drag; // drag coefficient
    float radius; // particle radius
    float m; // mass (kg)
    float invm; // 1 / mass (1/kg)
    float cr; // restitution
    float boxw; // simulation box width
    float boxh; // simulation box height
} pconfig_t;

typedef void (*pcallback_t)(void *ctx);

typedef struct {
    size_t count;
    size_t capacity;
    pcallback_t *items;
} pforcelist_t;

typedef struct { float vec[DIMS]; } pvec_t;

typedef struct {
    pconfig_t config;
    float time; // simulation runtime

    int count; // particle count

    pvec_t *x; // current state position
    pvec_t *v; // current state velocity
    pvec_t *xdot; // state position delta
    pvec_t *vdot; // state velocity delta
    pvec_t *workx; // position work vector
    pvec_t *workv; // velocity work vector
    pvec_t *f; // particle force

    pforcelist_t forcecallbacks;
} pstate_t;

void psys_init(pstate_t *psys, pconfig_t config, int count);
void psys_delete(pstate_t *psys);
void psys_step(pstate_t *psys, float delta_time);

#endif // PARTICLES_H
