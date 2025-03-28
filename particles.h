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

struct psysvec { float vec[DIMS]; };
struct pstatevec {
    struct psysvec x;
    struct psysvec v;
};

struct psys {
    struct psysconfig config;
    float time; // simulation runtime

    int count; // particle count

    struct psysvec *statex; // current state position
    struct psysvec *statev; // current state velocity
    struct psysvec *dstatex; // state position delta
    struct psysvec *dstatev; // state velocity delta
    struct psysvec *workx; // position work vector
    struct psysvec *workv; // velocity work vector
    struct psysvec *force; // particle force

    struct pforcelist forcecallbacks;
};

void psys_init(struct psys *psys, struct psysconfig config, int count);
void psys_delete(struct psys *psys);
void psys_step(struct psys *psys, float delta_time);

#endif // PARTICLES_H
