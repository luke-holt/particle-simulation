#ifndef PVEC_H
#define PVEC_H

#include <math.h>

typedef struct { float x, y; } pvec_t;

static inline void pvec_accum(pvec_t *accum, pvec_t vec) {
    accum->x += vec.x;
    accum->y += vec.y;
}

static inline void pvec_accum_many(int count, pvec_t *accum, const pvec_t *vec) {
    while (count--) pvec_accum(accum++, *(vec++));
}

static inline pvec_t pvec_scale(pvec_t vec, float d) {
    vec.x *= d;
    vec.y *= d;
    return vec;
}

static inline void pvec_scale_many(int count, pvec_t *out, pvec_t *in, float d) {
    while (count--) *(out++) = pvec_scale(*(in++), d);
}

static inline pvec_t pvec_sub(pvec_t a, pvec_t b) {
    return (pvec_t){ a.x-b.x, a.y-b.y };
}

static inline pvec_t pvec_add(pvec_t a, pvec_t b) {
    return (pvec_t){ a.x+b.x, a.y+b.y };
}

static inline float pvec_magsq(pvec_t v) {
    return v.x*v.x+v.y*v.y;
}

static inline float pvec_mag(pvec_t v) {
    return sqrtf(pvec_magsq(v));
}

static inline pvec_t pvec_normalize(pvec_t v) {
    return pvec_scale(v, 1.0 / pvec_mag(v));
}

static inline pvec_t pvec_dot(pvec_t a, pvec_t b) {
    return (pvec_t){ a.x * b.x, a.y * b.y };
}

static inline pvec_t pvec_negate(pvec_t v) {
    return pvec_scale(v, -1.0);
}

#endif // PVEC_H
