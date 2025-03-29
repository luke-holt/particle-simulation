#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>

#include "particles.h"

#include "util.h"

#ifndef PASSERT
#include <assert.h>
#define PASSERT assert
#endif

static inline void pvec_accum(pvec_t *accum, pvec_t vec) {
    accum->vec[0] += vec.vec[0];
    accum->vec[1] += vec.vec[1];
}
static inline void pvec_accum_many(int count, pvec_t *accum, const pvec_t *vec) {
    while (count--) pvec_accum(accum++, *(vec++));
}
static inline pvec_t pvec_scale(pvec_t vec, float d) {
    vec.vec[0] *= d;
    vec.vec[1] *= d;
    return vec;
}
static inline void pvec_scale_many(int count, pvec_t *out, pvec_t *in, float d) {
    while (count--) *(out++) = pvec_scale(*(in++), d);
}
static inline pvec_t pvec_rand(pvec_t range) {
    return (pvec_t) {
        ((float)rand() / (float)RAND_MAX) * range.vec[0],
        ((float)rand() / (float)RAND_MAX) * range.vec[1],
    };
}
static inline pvec_t pvec_sub(pvec_t a, pvec_t b) {
    return (pvec_t){ a.vec[0]-b.vec[0], a.vec[1]-b.vec[1] };
}
static inline pvec_t pvec_add(pvec_t a, pvec_t b) {
    return (pvec_t){ a.vec[0]+b.vec[0], a.vec[1]+b.vec[1] };
}
static inline float pvec_magsq(pvec_t v) {
    return v.vec[0]*v.vec[0]+v.vec[1]*v.vec[1];
}
static inline float pvec_mag(pvec_t v) {
    return sqrtf(pvec_magsq(v));
}
static inline pvec_t pvec_normalize(pvec_t v) {
    return pvec_scale(v, 1.0 / pvec_mag(v));
}
// particle collision
static inline bool pcol(float r, pvec_t xa, pvec_t xb) {
    return (r * r) > pvec_magsq(pvec_sub(xa, xb));
}

static void calculate_forces(pstate_t *psys);
static void derivative(pstate_t *psys);
static bool collisions(pstate_t *psys, float delta_time, int *p, int *q, float *coltime);
static void handle_collision(pstate_t *psys, float coltime, int p, int q);
static void step(pstate_t *psys, float delta_time);
static void wall_collisions(pstate_t *sys);

void
particles_init(pstate_t *psys, pconfig_t config, int count)
{
    PASSERT(psys);
    PASSERT(count);

    psys->count = count;
    psys->config = config;

    // alloc buffers
    psys->x = malloc(count*sizeof(*psys->x));
    PASSERT(psys->x);
    psys->v = malloc(count*sizeof(*psys->v));
    PASSERT(psys->v);
    psys->xdot = malloc(count*sizeof(*psys->xdot));
    PASSERT(psys->xdot);
    psys->vdot = malloc(count*sizeof(*psys->vdot));
    PASSERT(psys->vdot);
    psys->f = malloc(count*sizeof(*psys->f));
    PASSERT(psys->f);

    da_init(&psys->callbacks);
    psys->time = 0.0;

    // random positions
    for (size_t i = 0; i < psys->count; i++)
        psys->x[i] = pvec_rand((pvec_t){config.boxw, config.boxh});
}

void
particles_register_cb(pstate_t *psys, pcallback_t cb)
{
    da_append(&psys->callbacks, cb);
}

void
particles_delete(pstate_t *psys)
{
    PASSERT(psys);
    free(psys->x);
    free(psys->v);
    free(psys->xdot);
    free(psys->vdot);
    free(psys->f);
    da_delete(&psys->callbacks);
    memset(psys, 0, sizeof(*psys));
}

void
particles_step(pstate_t *psys, float delta_time)
{
    // clear forces
    memset(psys->f, 0, sizeof(*psys->f)*psys->count);

    calculate_forces(psys);

    wall_collisions(psys);

    derivative(psys);

    // apply step

    float coltime;
    int p, q;
    while (collisions(psys, delta_time, &p, &q, &coltime)) {

        // step until collision (x + xdot * coltime)
        step(psys, coltime);
        delta_time -= coltime;

        handle_collision(psys, coltime, p, q);

        derivative(psys);
    }

    // step the remainder
    step(psys, delta_time);
}

void
step(pstate_t *psys, float delta_time)
{
    for (int i = 0; i < psys->count; i++) {
        pvec_accum(&psys->x[i], pvec_scale(psys->xdot[i], delta_time));
        pvec_accum(&psys->v[i], pvec_scale(psys->vdot[i], delta_time));
    }
    psys->time += delta_time;
}

void
calculate_forces(pstate_t *psys)
{
    // gravity
    for (size_t i = 0; i < psys->count; i++) {
        psys->f[i].vec[1] += psys->config.gravity * psys->config.m;
    }

    // drag
    for (size_t i = 0; i < psys->count; i++)
        pvec_accum(&psys->f[i], pvec_scale(psys->v[i], psys->config.drag));

    // additional forces
    for (size_t i = 0; i < psys->callbacks.count; i++)
        psys->callbacks.items[i](psys);
}

void
derivative(pstate_t *psys)
{
    for (size_t i = 0; i < psys->count; i++) {
        // xdot = v
        psys->xdot[i] = psys->v[i];
        // vdot = f/m
        psys->vdot[i] = pvec_scale(psys->f[i], psys->config.invm);
    }
}

bool
collisions(pstate_t *psys, float delta_time, int *p, int *q, float *coltime)
{
    PASSERT(psys);
    PASSERT(p);
    PASSERT(q);
    PASSERT(coltime);

    *coltime = delta_time;
    bool c = false;
    for (int i = 0; i < psys->count; i++) {
        for (int j = 0; j < psys->count; j++) {
            if (i == j)
                continue;

            pvec_t xi, xj;
            // x + xdot * dt
            xi = pvec_add(psys->x[i], pvec_scale(psys->xdot[i], delta_time));
            xj = pvec_add(psys->x[j], pvec_scale(psys->xdot[j], delta_time));

            if (!pcol(psys->config.radius * 2, xi, xj))
                continue;

            c = true;

            float timestep = delta_time * 0.5;
            float offset = timestep * 0.5;

            int n = 10;
            do {
                // x + xdot * dt
                xi = pvec_add(psys->x[i], pvec_scale(psys->xdot[i], timestep));
                xj = pvec_add(psys->x[j], pvec_scale(psys->xdot[j], timestep));

                if (pcol(psys->config.radius * 2, xi, xj)) {
                    timestep -= offset;
                } else {
                    timestep += offset;
                }

                offset *= 0.5;
            } while (n--);

            if (*coltime > timestep) {
                *coltime = timestep;
                *p = i;
                *q = j;
            }
        }
    }

    return c;
}

void
handle_collision(pstate_t *psys, float coltime, int p, int q)
{
    // v2 = (1 + e)vcom - ev1
    // vcom = (vq + vq) / 2 <- for equal mass
    pvec_t vcom = pvec_scale(pvec_add(psys->v[p], psys->v[q]), 0.5 * (1 + psys->config.cr));
    pvec_t evp = pvec_scale(psys->v[p], psys->config.cr);
    pvec_t evq = pvec_scale(psys->v[q], psys->config.cr);
    psys->v[p] = pvec_sub(vcom, evp);
    psys->v[q] = pvec_sub(vcom, evq);

    // nudge particles apart if still colliding
    if (pcol(psys->config.radius*2, psys->x[q], psys->x[p])) {
        pvec_t dx = pvec_sub(psys->x[q], psys->x[p]);
        pvec_t dxu = pvec_normalize(dx);
        float delta = (psys->config.radius - pvec_mag(dx)) * 0.5;
        pvec_accum(&psys->x[p], pvec_scale(dxu, delta));
        pvec_accum(&psys->x[q], pvec_scale(dxu, -delta));
    }
}

void
wall_collisions(pstate_t *sys)
{
    for (size_t i = 0; i < sys->count; i++) {
        if ((sys->x[i].vec[0] - sys->config.radius) <= 0.0 && sys->v[i].vec[0] < 0.0 ||
            (sys->x[i].vec[0] + sys->config.radius) >= sys->config.boxw && sys->v[i].vec[0] > 0.0) {
            sys->v[i].vec[0] *= -1.0 * sys->config.cr;
        }
        if ((sys->x[i].vec[1] - sys->config.radius) <= 0.0 && sys->v[i].vec[1] < 0.0 ||
            (sys->x[i].vec[1] + sys->config.radius) >= sys->config.boxh && sys->v[i].vec[1] > 0.0) {
            sys->v[i].vec[1] *= -1.0 * sys->config.cr;
        }
    }
}

