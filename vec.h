#include <math.h>

typedef struct { float x, y; } vec2;

#define as_vec2(f) (*(vec2 *)((void *)&(f)))

static inline vec2 sub(vec2 a, vec2 b) { return (vec2){ a.x-b.x, a.y-b.y }; }
static inline vec2 add(vec2 a, vec2 b) { return (vec2){ a.x+b.x, a.y+b.y }; }
static inline float dot(vec2 a, vec2 b) { return a.x*b.x+a.y*b.y; }
static inline float magsq(vec2 v) { return v.x*v.x+v.y*v.y; }
static inline float mag(vec2 v) { return sqrtf(magsq(v)); }
static inline vec2 scale(vec2 v, float s) { return (vec2){ v.x*s, v.y*s }; }
static inline vec2 cast(vec2 a, vec2 b) { return scale(b, dot(a, b) / magsq(b)); }

