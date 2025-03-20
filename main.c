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

#define GRAVITY (-9.81)

const int fps = 60;
const double delta_time = 1.0 / fps;
const int scw = 800;
const int sch = 600;

float magnitude(float v[2]) { return sqrtf(v[0]*v[0]+v[1]*v[1]); }

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

void spring(void *ctx) {
    ParticleSystem *sys = (ParticleSystem *)ctx;
    for (size_t i = 0; i < sys->particles.count; i++) {
            for (size_t j = 0; j < sys->particles.count; j++) {
            if (i == j) continue;
            Particle *a = &sys->particles.items[i];
            Particle *b = &sys->particles.items[j];

            float dx[2], dv[2], mdx, f[2];
            dx[0] = a->x[0] - b->x[0];
            dx[1] = a->x[1] - b->x[1];
            dv[0] = a->v[0] - b->v[0];
            dv[1] = a->v[1] - b->v[1];
            mdx = magnitude(dx);
            float ks = 1;
            float kd = 1;

            f[0] = -(ks * (mdx - 30) + kd * dv[0] * dx[0] / mdx) * dx[0] / mdx;
            f[1] = -(ks * (mdx - 30) + kd * dv[1] * dx[1] / mdx) * dx[1] / mdx;

            a->f[0] += f[0];
            a->f[1] += f[1];
            b->f[0] -= f[0];
            b->f[1] -= f[1];
        }
    }
}


int
main(void)
{
    ParticleSystem sys = particle_system_new(10, scw, sch, 10, 0.0, 0.10);

    // da_append(&sys.forces, collision);
    da_append(&sys.forces, spring);

    InitWindow(scw, sch, "physical-modelling-particles");
    SetTargetFPS(60);

    while (!WindowShouldClose()) {

        euler_step(&sys, delta_time);

        BeginDrawing();
        ClearBackground(BLACK);

        for (size_t i = 0; i < sys.particles.count; i++) {
            for (size_t j = 0; j < sys.particles.count; j++) {
                if (i == j) continue;
                Particle *a = &sys.particles.items[i];
                Particle *b = &sys.particles.items[j];
                Color c = (Color) {20, 20, 20, 255};
                DrawLineV(*(Vector2 *)a->x, *(Vector2 *)b->x, c);
            }
        }
        for (size_t i = 0; i < sys.particles.count; i++) {
            // DrawPixelV(*(Vector2 *)sys.particles.items[i].x, RAYWHITE);
            DrawCircleV(*(Vector2 *)sys.particles.items[i].x, 3.0, RAYWHITE);
        }

        EndDrawing();
    }

    return 0;
}

