#ifndef PVEC_H
#define PVEC_H

#include <math.h>

#define DIMS (2)

typedef struct { float vec[DIMS]; } pvec_t;

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

static inline pvec_t pvec_dot(pvec_t a, pvec_t b) {
    return (pvec_t){ a.vec[0] * b.vec[0], a.vec[1] * b.vec[1] };
}

static inline pvec_t pvec_negate(pvec_t v) {
    return pvec_scale(v, -1.0);
}

#endif // PVEC_H
