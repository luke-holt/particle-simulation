#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>

#include "particles.h"

#include "util.h"

#ifndef P2ASSERT
#include <assert.h>
#define P2ASSERT assert
#endif

static void pvec_accum(struct pvec *accum, struct pvec vec) {
    accum->vec[0] += vec.vec[0];
    accum->vec[1] += vec.vec[1];
}
static void pvec_accum_many(int count, struct pvec *accum, const struct pvec *vec) {
    while (count--) pvec_accum(accum++, *(vec++));
}
static struct pvec pvec_scale(struct pvec vec, float d) {
    vec.vec[0] *= d;
    vec.vec[1] *= d;
    return vec;
}
static void pvec_scale_many(int count, struct pvec *out, struct pvec *in, float d) {
    while (count--) *(out++) = pvec_scale(*(in++), d);
}
static struct pvec pvec_rand(struct pvec range) {
    return (struct pvec) {
        ((float)rand() / (float)RAND_MAX) * range.vec[0],
        ((float)rand() / (float)RAND_MAX) * range.vec[1],
    };
}
static struct pvec pvec_sub(struct pvec a, struct pvec b) {
    return (struct pvec){ a.vec[0]-b.vec[0], a.vec[1]-b.vec[1] };
}
static struct pvec pvec_add(struct pvec a, struct pvec b) {
    return (struct pvec){ a.vec[0]+b.vec[0], a.vec[1]+b.vec[1] };
}
static inline float pvec_magsq(struct pvec v) {
    return v.vec[0]*v.vec[0]+v.vec[1]*v.vec[1];
}
static inline float pvec_mag(struct pvec v) {
    return sqrtf(pvec_magsq(v));
}
static inline struct pvec pvec_normalize(struct pvec v) {
    return pvec_scale(v, 1.0 / pvec_mag(v));
}

// particle collision
static bool pcol(float r, struct pvec xa, struct pvec xb) {
    return (r * r) > pvec_magsq(pvec_sub(xa, xb));
}

static void calculate_forces2(struct psys *psys);
static void derivative(struct psys *psys);
bool collisions(struct psys *psys, float delta_time, int *p, int *q, float *coltime);
void handle_collision(struct psys *psys, float coltime, int p, int q);

void
psys_init(struct psys *psys, struct psysconfig config, int count)
{
    P2ASSERT(psys);
    P2ASSERT(count);

    psys->count = count;
    psys->config = config;

    // alloc buffers
    psys->x = malloc(count*sizeof(*psys->x));
    P2ASSERT(psys->x);
    psys->v = malloc(count*sizeof(*psys->v));
    P2ASSERT(psys->v);
    psys->workx = malloc(count*sizeof(*psys->workx));
    P2ASSERT(psys->workx);
    psys->workv = malloc(count*sizeof(*psys->workv));
    P2ASSERT(psys->workv);
    psys->xdot = malloc(count*sizeof(*psys->xdot));
    P2ASSERT(psys->xdot);
    psys->vdot = malloc(count*sizeof(*psys->vdot));
    P2ASSERT(psys->vdot);
    psys->f = malloc(count*sizeof(*psys->f));
    P2ASSERT(psys->f);

    da_init(&psys->forcecallbacks);
    psys->time = 0.0;

    // random positions
    for (size_t i = 0; i < psys->count; i++)
        psys->x[i] = pvec_rand((struct pvec){config.boxw, config.boxh});
}

void
psys_delete(struct psys *psys)
{
    P2ASSERT(psys);
    free(psys->x);
    free(psys->v);
    free(psys->xdot);
    free(psys->vdot);
    free(psys->workx);
    free(psys->workv);
    free(psys->f);
    da_delete(&psys->forcecallbacks);
    memset(psys, 0, sizeof(*psys));
}

void
psys_step(struct psys *psys, float delta_time)
{
    // clear forces
    memset(psys->f, 0, sizeof(*psys->f)*psys->count);

    calculate_forces2(psys);

    derivative(psys);

    // apply step

    float duration = delta_time;
    float coltime;
    int p, q;
    while (collisions(psys, duration, &p, &q, &coltime)) {

        // step until collision (x + xdot * coltime)
        pvec_scale_many(psys->count, psys->workx, psys->xdot, coltime);
        pvec_scale_many(psys->count, psys->workv, psys->vdot, coltime);
        pvec_accum_many(psys->count, psys->x, psys->workx);
        pvec_accum_many(psys->count, psys->v, psys->workv);
        duration -= coltime;

        handle_collision(psys, coltime, p, q);

        derivative(psys);
    }

    // step the remainder
    pvec_scale_many(psys->count, psys->workx, psys->xdot, duration);
    pvec_scale_many(psys->count, psys->workv, psys->vdot, duration);
    pvec_accum_many(psys->count, psys->x, psys->workx);
    pvec_accum_many(psys->count, psys->v, psys->workv);

    // scale derivative by timestep

    // pvec_scale_many(psys->count, psys->workx, psys->xdot, delta_time);
    // pvec_scale_many(psys->count, psys->workv, psys->vdot, delta_time);

    // pvec_accum_many(psys->count, psys->x, psys->workx);
    // pvec_accum_many(psys->count, psys->v, psys->workv);

    psys->time += delta_time;
}

void
calculate_forces2(struct psys *psys)
{
    // gravity
    for (size_t i = 0; i < psys->count; i++) {
        psys->f[i].vec[1] += psys->config.gravity * psys->config.m;
    }

    // drag
    for (size_t i = 0; i < psys->count; i++)
        pvec_accum(&psys->f[i], pvec_scale(psys->v[i], psys->config.drag));

    // additional forces
    for (size_t i = 0; i < psys->forcecallbacks.count; i++)
        psys->forcecallbacks.items[i](psys);
}

void
derivative(struct psys *psys)
{
    for (size_t i = 0; i < psys->count; i++) {
        // xdot = v
        psys->xdot[i] = psys->v[i];
        // vdot = f/m
        psys->vdot[i] = pvec_scale(psys->f[i], psys->config.invm);
    }
}

bool
collisions(struct psys *psys, float delta_time, int *p, int *q, float *coltime)
{
    P2ASSERT(psys);
    P2ASSERT(p);
    P2ASSERT(q);
    P2ASSERT(coltime);

    *coltime = delta_time;
    bool c = false;
    for (int i = 0; i < psys->count; i++) {
        for (int j = 0; j < psys->count; j++) {
            if (i == j)
                continue;

            struct pvec xi, xj;
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
handle_collision(struct psys *psys, float coltime, int p, int q)
{
    // v2 = (1 + e)vcom - ev1
    // vcom = (vq + vq) / 2 <- for equal mass
    struct pvec vcom = pvec_scale(pvec_add(psys->v[p], psys->v[q]), 0.5 * (1 + psys->config.cr));
    struct pvec evp = pvec_scale(psys->v[p], psys->config.cr);
    struct pvec evq = pvec_scale(psys->v[q], psys->config.cr);
    psys->v[p] = pvec_sub(vcom, evp);
    psys->v[q] = pvec_sub(vcom, evq);

    // nudge particles apart if still colliding
    if (pcol(psys->config.radius*2, psys->x[q], psys->x[p])) {
        struct pvec dx = pvec_sub(psys->x[q], psys->x[p]);
        struct pvec dxu = pvec_normalize(dx);
        float delta = (psys->config.radius - pvec_mag(dx)) * 0.5;
        pvec_accum(&psys->x[p], pvec_scale(dxu, delta));
        pvec_accum(&psys->x[q], pvec_scale(dxu, -delta));
    }
}
