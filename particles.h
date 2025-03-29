#ifndef PARTICLES_H
#define PARTICLES_H

#include <stddef.h>

#define DIMS (2)

struct psysconfig {
    float gravity; // system gravity (m/s^2)
    float drag; // drag coefficient
    float radius; // particle radius
    float m; // mass (kg)
    float invm; // 1 / mass (1/kg)
    float cr; // restitution
    float boxw; // simulation box width
    float boxh; // simulation box height
};

typedef void (*psys_force_callback)(void *ctx);

struct pforcelist {
    size_t count;
    size_t capacity;
    psys_force_callback *items;
};

struct pvec { float vec[DIMS]; };
struct pstatevec {
    struct pvec x;
    struct pvec v;
};

struct psys {
    struct psysconfig config;
    float time; // simulation runtime

    int count; // particle count

    struct pvec *x; // current state position
    struct pvec *v; // current state velocity
    struct pvec *xdot; // state position delta
    struct pvec *vdot; // state velocity delta
    struct pvec *workx; // position work vector
    struct pvec *workv; // velocity work vector
    struct pvec *f; // particle force

    struct pforcelist forcecallbacks;
};

void psys_init(struct psys *psys, struct psysconfig config, int count);
void psys_delete(struct psys *psys);
void psys_step(struct psys *psys, float delta_time);

#endif // PARTICLES_H
