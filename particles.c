#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include "particles.h"

#include "util.h"

#ifndef P2ASSERT
#include <assert.h>
#define P2ASSERT assert
#endif

static void psysvec_accum(struct psysvec *accum, struct psysvec vec) {
    accum->vec[0] += vec.vec[0];
    accum->vec[1] += vec.vec[1];
}
static void psysvec_accum_many(int count, struct psysvec *accum, const struct psysvec *vec) {
    while (count--) psysvec_accum(accum++, *(vec++));
}
static struct psysvec psysvec_scale(struct psysvec vec, float d) {
    vec.vec[0] *= d;
    vec.vec[1] *= d;
    return vec;
}
static void psysvec_scale_many(int count, struct psysvec *out, struct psysvec *in, float d) {
    while (count--) *(out++) = psysvec_scale(*(in++), d);
}
static struct psysvec psysvec_rand(struct psysvec range) {
    return (struct psysvec) {
        ((float)rand() / (float)RAND_MAX) * range.vec[0],
        ((float)rand() / (float)RAND_MAX) * range.vec[1],
    };
}
static struct psysvec psysvec_sub(struct psysvec a, struct psysvec b) {
    return (struct psysvec){ a.vec[0]-b.vec[0], a.vec[1]-b.vec[1] };
}
static struct psysvec psysvec_add(struct psysvec a, struct psysvec b) {
    return (struct psysvec){ a.vec[0]+b.vec[0], a.vec[1]+b.vec[1] };
}
static inline float psysvec_magsq(struct psysvec v) {
    return v.vec[0]*v.vec[0]+v.vec[1]*v.vec[1];
}
static inline float psysvec_mag(struct psysvec v) {
    return sqrtf(psysvec_magsq(v));
}

// particle collision
static bool pcol(float radius, struct psysvec xa, struct psysvec xb) {
    return (radius * radius) < psysvec_magsq(psysvec_sub(xa, xb));
}

static void calculate_forces2(struct psys *psys);
static void derivative(struct psys *psys);

void
psys_init(struct psys *psys, struct psysconfig config, int count)
{
    P2ASSERT(psys);
    P2ASSERT(count);

    psys->count = count;
    psys->config = config;

    // alloc buffers
    psys->statex = malloc(count*sizeof(*psys->statex));
    P2ASSERT(psys->statex);
    psys->statev = malloc(count*sizeof(*psys->statev));
    P2ASSERT(psys->statev);
    psys->workx = malloc(count*sizeof(*psys->workx));
    P2ASSERT(psys->workx);
    psys->workv = malloc(count*sizeof(*psys->workv));
    P2ASSERT(psys->workv);
    psys->dstatex = malloc(count*sizeof(*psys->dstatex));
    P2ASSERT(psys->dstatex);
    psys->dstatev = malloc(count*sizeof(*psys->dstatev));
    P2ASSERT(psys->dstatev);
    psys->force = malloc(count*sizeof(*psys->force));
    P2ASSERT(psys->force);

    da_init(&psys->forcecallbacks);
    psys->time = 0.0;

    // random positions
    for (size_t i = 0; i < psys->count; i++)
        psys->statex[i] = psysvec_rand((struct psysvec){config.boxw, config.boxh});
}

void
psys_delete(struct psys *psys)
{
    P2ASSERT(psys);
    free(psys->statex);
    free(psys->statev);
    free(psys->workx);
    free(psys->workv);
    free(psys->dstatex);
    free(psys->dstatev);
    free(psys->force);
    da_delete(&psys->forcecallbacks);
    memset(psys, 0, sizeof(*psys));
}

void
psys_step(struct psys *psys, float delta_time)
{
    // clear forces
    memset(psys->force, 0, sizeof(*psys->force)*psys->count);

    calculate_forces2(psys);

    derivative(psys);

    // apply step

    // scale derivative by timestep

    psysvec_scale_many(psys->count, psys->workx, psys->dstatex, delta_time);
    psysvec_scale_many(psys->count, psys->workv, psys->dstatev, delta_time);

#if 0
    for (int i = 0; i < psys->count; i++) {
        for (int j = 0; j < psys->count; j++) {
            if (i == j)
                continue;

            struct psysvec xi = psysvec_add(psys->statex[i], psysvec_scale(psys->dstatex[i], delta_time));
            struct psysvec xj = psysvec_add(psys->statex[j], psysvec_scale(psys->dstatex[j], delta_time));
            struct psysvec vi = psysvec_add(psys->statev[i], psysvec_scale(psys->dstatev[i], delta_time));
            struct psysvec vj = psysvec_add(psys->statev[j], psysvec_scale(psys->dstatev[j], delta_time));

            // particle collision
            if (!pcol(psys->config.radius, xi, xj))
                continue;

            // collision time?
            struct psysvec wx = psysvec_scale(psys->dstatex[i], delta_time);
            struct psysvec wv = psysvec_scale(psys->dstatev[i], delta_time);


        }
    }
#endif

    psysvec_accum_many(psys->count, psys->statex, psys->workx);
    psysvec_accum_many(psys->count, psys->statev, psys->workv);

    psys->time += delta_time;
}

void
calculate_forces2(struct psys *psys)
{
    // gravity
    for (size_t i = 0; i < psys->count; i++) {
        psys->force[i].vec[1] += psys->config.gravity * psys->config.m;
    }

    // drag
    for (size_t i = 0; i < psys->count; i++)
        psysvec_accum(&psys->force[i], psysvec_scale(psys->statev[i], psys->config.drag));

    // additional forces
    for (size_t i = 0; i < psys->forcecallbacks.count; i++)
        psys->forcecallbacks.items[i](psys);
}

void
derivative(struct psys *psys)
{
    for (size_t i = 0; i < psys->count; i++) {
        // xdot = v
        psys->dstatex[i] = psys->statev[i];
        // vdot = f/m
        psys->dstatev[i] = psysvec_scale(psys->force[i], psys->config.invm);
    }
}
