#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "raylib.h"

#define PARTICLES_IMPLEMENTATION
#include "particles.h"

#define UTIL_IMPLEMENTATION
#include "util.h"

/*

f = ma
xdot2 = f/m <- second derivative of position
to reduce to first order, introduce velocity v = xdot1

vdot1 = f/m, xdot1 = v

phase space equation of motion (2 dimensions)

[xdot1_1, xdot1_2, vdot1_1, vdot1_2]
=
[v1, v2, f1/m, f2/m]

*/

const int fps = 60;
const double delta_time = 1.0 / fps;
const int scw = 800;
const int sch = 600;

void gravity(void *ctx) {
    ParticleSystem *sys = (ParticleSystem *)ctx;
    for (size_t i = 0; i < sys->particles.count; i++) {
        double m = sys->particles.items[i].m;
        sys->particles.items[i].f[1] += -9.81 * m;
    }
}

void drag(void *ctx) {
    ParticleSystem *sys = (ParticleSystem *)ctx;
    for (size_t i = 0; i < sys->particles.count; i++) {
        double d[2];
        d[0] = sys->particles.items[i].v[0] * -0.10;
        d[1] = sys->particles.items[i].v[1] * -0.10;
        sys->particles.items[i].f[0] += d[0];
        sys->particles.items[i].f[1] += d[1];
    }
}

void collision(void *ctx) {
    ParticleSystem *sys = (ParticleSystem *)ctx;
    for (size_t i = 0; i < sys->particles.count; i++) {
        if (sys->particles.items[i].x[0] <= 0.0)
            sys->particles.items[i].v[0] *= -1.0;
        if (sys->particles.items[i].x[1] <= 0.0)
            sys->particles.items[i].v[1] *= -1.0;
        if (sys->particles.items[i].x[0] >= scw)
            sys->particles.items[i].v[0] *= -1.0;
        if (sys->particles.items[i].x[1] >= sch)
            sys->particles.items[i].v[1] *= -1.0;
    }
}

float magnitude(float v[2]) {
    return sqrtf(v[0]*v[0]+v[1]*v[1]);
}

int
main(void)
{
    ParticleSystem sys = particle_system_new(10000, scw, sch, 10);

    da_append(&sys.forces, gravity);
    da_append(&sys.forces, drag);
    da_append(&sys.forces, collision);

    InitWindow(scw, sch, "physical-modelling-particles");
    SetTargetFPS(60);

    while (!WindowShouldClose()) {

        euler_step(&sys, delta_time);

        BeginDrawing();
        ClearBackground(BLACK);

        for (size_t i = 0; i < sys.particles.count; i++) {
            DrawPixelV(*(Vector2 *)sys.particles.items[i].x, BLUE);
        }

        EndDrawing();
    }

    return 0;
}

