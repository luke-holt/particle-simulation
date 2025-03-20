#ifndef PARTICLES_H
#define PARTICLES_H

#include <stddef.h>

// dimensions
#define DIMS (2)

typedef struct {
    float m;       // mass
    float x[DIMS]; // position vector
    float v[DIMS]; // velocity vector
    float f[DIMS]; // force accumulator
} Particle;

typedef struct {
    size_t count;
    size_t capacity;
    Particle *items;
} ParticleList;

typedef void (*force_callback)(void *ctx);

typedef struct {
    size_t count;
    size_t capacity;
    force_callback *items;
} ForceList;

typedef struct {
    float t;    // simulation clock
    float gravity;
    float drag;
    float *dbuf;
    float *statebuf;
    ParticleList particles;
    ForceList forces;
} ParticleSystem;

ParticleSystem particle_system_new(int count, int width, int height, float mass, float gravity, float drag);
void particle_system_delete(ParticleSystem *sys);

#endif // PARTICLES_H

#define PARTICLES_IMPLEMENTATION
#ifdef PARTICLES_IMPLEMENTATION

#include <stdlib.h>
#include <string.h>

#ifndef PARTICLES_ASSERT
#include <assert.h>
#define PARTICLES_ASSERT assert
#endif // PARTICLES_ASSERT

#include "util.h"

void calculate_forces(ParticleSystem *sys);
int particle_dims(ParticleSystem *sys);
void particle_get_state(ParticleSystem *sys, float *dst);
void particle_set_state(ParticleSystem *sys, float *src);
void clear_forces(ParticleSystem *sys);
void particle_derivative(ParticleSystem *sys, float *dst);
void scale_vector(float *vector, int n, float d);
void add_vector(float *v1, float *v2, float *out, int n);
void euler_step(ParticleSystem *sys, float delta_time);

// length of state derivative and force vectors
int
particle_dims(ParticleSystem *sys)
{
    // 2 n-dimensional state vectors for each particle
    return (DIMS * 2 * sys->particles.count);
}

// gather state from the particles into dst
void
particle_get_state(ParticleSystem *sys, float *dst)
{
    for (int i = 0; i < sys->particles.count; i++) {
        *(dst++) = sys->particles.items[i].x[0];
        *(dst++) = sys->particles.items[i].x[1];
        *(dst++) = sys->particles.items[i].v[0];
        *(dst++) = sys->particles.items[i].v[1];
    }
}

// scatter state from src into the particles
void
particle_set_state(ParticleSystem *sys, float *src)
{
    for (int i = 0; i < sys->particles.count; i++) {
        sys->particles.items[i].x[0] = *(src++);
        sys->particles.items[i].x[1] = *(src++);
        sys->particles.items[i].v[0] = *(src++);
        sys->particles.items[i].v[1] = *(src++);
    }
}

void
clear_forces(ParticleSystem *sys)
{
    for (int i = 0; i < sys->particles.count; i++) {
        sys->particles.items[i].f[0] = 0.0;
        sys->particles.items[i].f[1] = 0.0;
    }
}

// calculate derivative, place in dst
void
particle_derivative(ParticleSystem *sys, float *dst)
{
    clear_forces(sys);
    calculate_forces(sys);
    for (int i = 0; i < sys->particles.count; i++) {
        // xdot = v
        *(dst++) = sys->particles.items[i].v[0];
        *(dst++) = sys->particles.items[i].v[1];
        // vdot = f/m
        *(dst++) = sys->particles.items[i].f[0] / sys->particles.items[i].m;
        *(dst++) = sys->particles.items[i].f[1] / sys->particles.items[i].m;
    }
}

void
scale_vector(float *vector, int n, float d)
{
    while (n--) *(vector++) *= d;
}

void
add_vector(float *v1, float *v2, float *out, int n)
{
    while (n--) *(out++) = *(v1++) + *(v2++);
}

void
euler_step(ParticleSystem *sys, float delta_time)
{
    int n = particle_dims(sys);
    particle_derivative(sys, sys->dbuf);
    scale_vector(sys->dbuf, n, delta_time);
    particle_get_state(sys, sys->statebuf);
    add_vector(sys->dbuf, sys->statebuf, sys->statebuf, n);
    particle_set_state(sys, sys->statebuf);
    sys->t += delta_time;
}

void
calculate_forces(ParticleSystem *sys)
{
    // gravity
    for (size_t i = 0; i < sys->particles.count; i++) {
        double m = sys->particles.items[i].m;
        sys->particles.items[i].f[1] += sys->gravity * m;
    }

    // drag
    for (size_t i = 0; i < sys->particles.count; i++) {
        sys->particles.items[i].f[0] += sys->particles.items[i].v[0] * sys->drag;
        sys->particles.items[i].f[1] += sys->particles.items[i].v[1] * sys->drag;
    }

    // additional forces
    for (size_t i = 0; i < sys->forces.count; i++)
        sys->forces.items[i](sys);
}

ParticleSystem
particle_system_new(int count, int width, int height, float mass, float gravity, float drag)
{
    ParticleSystem sys;
    sys.t = 0.0;
    sys.gravity = gravity;
    sys.drag = drag;
    da_init(&sys.forces);
    da_init(&sys.particles);
    da_resize(&sys.particles, count);
    sys.particles.count = count;
    memset(sys.particles.items, 0, sizeof(*sys.particles.items)*count);
    // scratch buffers used in euler_step
    sys.dbuf = (float *)malloc(particle_dims(&sys) * sizeof(*sys.dbuf));
    PARTICLES_ASSERT(sys.dbuf);
    sys.statebuf = (float *)malloc(particle_dims(&sys) * sizeof(*sys.statebuf));
    PARTICLES_ASSERT(sys.statebuf);
    // random positions
    for (size_t i = 0; i < sys.particles.count; i++) {
        sys.particles.items[i].x[0] = ((float)rand() / (float)RAND_MAX) * width;
        sys.particles.items[i].x[1] = ((float)rand() / (float)RAND_MAX) * height;
        sys.particles.items[i].m = mass;
    }
    return sys;
}

void
particle_system_delete(ParticleSystem *sys)
{
    assert(sys);
    da_delete(&sys->particles);
    da_delete(&sys->forces);
    free(sys->dbuf);
    free(sys->statebuf);
    memset(sys, 0, sizeof(*sys));
}

#endif // PARTICLES_IMPLEMENTATION
