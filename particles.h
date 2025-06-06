#ifndef PARTICLES_H
#define PARTICLES_H

#include <stddef.h>

#include "pvec.h"

typedef struct {
    pvec_t gravity; // system gravity (m/s^2)
    pvec_t box; // simulation box
    float drag; // drag coefficient
    float radius; // particle radius
    float m; // mass (kg)
    float invm; // 1 / mass (1/kg)
    float cr; // restitution
} pconfig_t;

typedef struct pstate_t pstate_t;

typedef void (*pcallback_t)(pstate_t *state);

typedef struct {
    size_t count;
    size_t capacity;
    pcallback_t *items;
} pcallbacklist_t;

struct pstate_t {
    pconfig_t config;
    float time; // simulation runtime

    int count; // particle count
    int capacity; // system capacity

    pvec_t *x; // current state position
    pvec_t *v; // current state velocity
    pvec_t *xdot; // state position delta
    pvec_t *vdot; // state velocity delta
    pvec_t *f; // particle force

    pcallbacklist_t callbacks;
};

void particles_init(pstate_t *psys, pconfig_t config, int capacity);
void particles_delete(pstate_t *psys);
void particles_step(pstate_t *psys, float delta_time);
void particles_register_cb(pstate_t *psys, pcallback_t cb);
void particles_add(pstate_t *psys, pvec_t x, pvec_t v);

#endif // PARTICLES_H
