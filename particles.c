#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>

#include "particles.h"

#include "util.h"

#ifndef PASSERT
#include <assert.h>
#define PASSERT assert
#endif

static inline pvec_t pvec_rand(pvec_t range) {
    return (pvec_t) {
        ((float)rand() / (float)RAND_MAX) * range.x,
        ((float)rand() / (float)RAND_MAX) * range.y,
    };
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

void
particles_init(pstate_t *psys, pconfig_t config, int capacity)
{
    PASSERT(psys);
    PASSERT(capacity);

    psys->capacity = capacity;
    psys->count = 0;
    psys->config = config;

    // alloc buffers
    psys->x = malloc(capacity*sizeof(*psys->x));
    PASSERT(psys->x);
    psys->v = malloc(capacity*sizeof(*psys->v));
    PASSERT(psys->v);
    psys->xdot = malloc(capacity*sizeof(*psys->xdot));
    PASSERT(psys->xdot);
    psys->vdot = malloc(capacity*sizeof(*psys->vdot));
    PASSERT(psys->vdot);
    psys->f = malloc(capacity*sizeof(*psys->f));
    PASSERT(psys->f);

    da_init(&psys->callbacks);
    psys->time = 0.0;
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
    memset(psys->f, 0, sizeof(*psys->f)*psys->capacity);

    calculate_forces(psys);

    // wall_collisions(psys);

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
    // gravity f = m*a
    for (size_t i = 0; i < psys->count; i++)
        pvec_accum(&psys->f[i], pvec_scale(psys->config.gravity, psys->config.m));

    // drag f_d = -v*d
    for (size_t i = 0; i < psys->count; i++)
        pvec_accum(&psys->f[i], pvec_scale(psys->v[i], -psys->config.drag));

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
    // assuming equal mass...
    // v1' = v1 - ((v1-v2) . (x1-x2)) / mag(x1-x2)^2 * (x1-x2)

    pvec_t dx = pvec_sub(psys->x[p], psys->x[q]);
    pvec_t dv = pvec_sub(psys->v[p], psys->v[q]);
    float dsq = pvec_magsq(dx);

    pvec_t vec = pvec_scale(pvec_dot(pvec_dot(dx, dv), dx), 1.0 / dsq);

    psys->v[p] = pvec_sub(psys->v[p], vec);
    psys->v[q] = pvec_sub(psys->v[q], pvec_negate(vec));

    // nudge particles apart if still colliding
    if (pcol(psys->config.radius*2, psys->x[q], psys->x[p])) {
        pvec_t dxu = pvec_normalize(dx);
        float delta = psys->config.radius - pvec_mag(dx);
        pvec_accum(&psys->x[p], pvec_scale(dxu, -delta*0.25));
        pvec_accum(&psys->x[q], pvec_scale(dxu, delta*0.25));
    }
}

void
wall_collisions(pstate_t *sys)
{
    for (size_t i = 0; i < sys->count; i++) {
        pvec_t *x = &sys->x[i];
        pvec_t *v = &sys->v[i];

        // left
        if (x->x - sys->config.radius < 0.0) {
            x->x = sys->config.radius;
            if (v->x < 0.0)
                v->x *= -sys->config.cr;
        }

        // right
        if (x->x + sys->config.radius > sys->config.box.x) {
            x->x = sys->config.box.x - sys->config.radius;
            if (v->x > 0.0)
                v->x *= -sys->config.cr;
        }

        // down
        if (x->y - sys->config.radius < 0.0) {
            x->y = sys->config.radius;
            if (v->y < 0.0)
                v->y *= -sys->config.cr;
        }

        // up
        if (x->y + sys->config.radius > sys->config.box.y) {
            x->y = sys->config.box.y - sys->config.radius;
            if (v->y > 0.0)
                v->y *= -sys->config.cr;
        }
    }
}

void
particles_add(pstate_t *psys, pvec_t x, pvec_t v)
{
    if (psys->count >= psys->capacity)
        return;

    psys->x[psys->count] = x;
    psys->v[psys->count] = v;
    psys->xdot[psys->count] = (pvec_t){0};
    psys->vdot[psys->count] = (pvec_t){0};
    psys->f[psys->count] = (pvec_t){0};
    psys->count++;
}
