#include <stdio.h>
#include <assert.h>

#include "raylib.h"

#include "particles.h"
#include "vec.h"

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

#define GRAVITY 0 // (-9.81)
#define DRAG (0.10)

const int fps = 60;
const float delta_time = 1.0 / fps;
const int scw = 800;
const int sch = 600;

#if 0
void wall_collisions(void *ctx) {
    ParticleSystem *sys = (ParticleSystem *)ctx;
    for (size_t i = 0; i < sys->particles.count; i++) {
        Particle *p = &sys->particles.items[i];
        if ((p->x[0] - sys->config.radius) <= 0.0 && p->v[0] < 0.0 ||
            (p->x[0] + sys->config.radius) >= scw && p->v[0] > 0.0) {
            p->v[0] *= -1.0 * sys->config.cr;
        }
        if ((p->x[1] - sys->config.radius) <= 0.0 && p->v[1] < 0.0 ||
            (p->x[1] + sys->config.radius) >= sch && p->v[1] > 0.0) {
            p->v[1] *= -1.0 * sys->config.cr;
        }
    }
}

void collisions(void *ctx) {
    ParticleSystem *sys = (ParticleSystem *)ctx;
    for (size_t i = 0; i < sys->particles.count; i++) {
        for (size_t j = 0; j < sys->particles.count; j++) {
            Particle *a = &sys->particles.items[i];
            Particle *b = &sys->particles.items[j];

            vec2 ax = as_vec2(a->x);
            vec2 bx = as_vec2(b->x);
            vec2 av = as_vec2(a->v);
            vec2 bv = as_vec2(b->v);

            // not colliding
            if ((i == j) || mag(sub(ax, bx)) > sys->config.radius)
                continue;

            float e = sys->config.cr;

            // velocity of center of mass (equal particle mass)
            vec2 vcom = scale(add(av, bv), 0.5);
            vec2 x = scale(vcom, (1 + e));

            av = sub(x, scale(av, e));
            bv = sub(x, scale(bv, e));

            a->v[0] = av.x; a->v[1] = av.y;
            b->v[0] = bv.x; b->v[1] = bv.y;
        }
    }
}
#endif

void spring(void *ctx) {
    struct psys *sys = (struct psys *)ctx;
    for (size_t i = 0; i < sys->count; i++) {
        for (size_t j = 0; j < sys->count; j++) {
            if (i == j) continue;

            float dx[2], dv[2], dxmag, f[2];
            dx[0] = sys->statex[i].vec[0] - sys->statex[j].vec[0];
            dx[1] = sys->statex[i].vec[1] - sys->statex[j].vec[1];
            dv[0] = sys->statev[i].vec[0] - sys->statev[j].vec[0];
            dv[1] = sys->statev[i].vec[1] - sys->statev[j].vec[1];

            dxmag = mag(as_vec2(dx));
            float ks = 1;
            float kd = 1;

            f[0] = -(ks * (dxmag - 30) + kd * dv[0] * dx[0] / dxmag) * dx[0] / dxmag;
            f[1] = -(ks * (dxmag - 30) + kd * dv[1] * dx[1] / dxmag) * dx[1] / dxmag;

            sys->force[i].vec[0] += f[0];
            sys->force[i].vec[1] += f[1];
            sys->force[j].vec[0] -= f[0];
            sys->force[j].vec[1] -= f[1];
        }
    }
}

void mouse_coupling(void *ctx) {
    struct psys *sys = (struct psys *)ctx;

    if (!IsMouseButtonDown(MOUSE_BUTTON_LEFT))
        return;

    Vector2 mouse_x = GetMousePosition();
    Vector2 mouse_v = GetMouseDelta();

    float dx[2], dv[2], dxmag, f[2];
    dx[0] = sys->statex[0].vec[0] - mouse_x.x;
    dx[1] = sys->statex[0].vec[1] - mouse_x.y;
    dv[0] = sys->statev[0].vec[0] - mouse_v.x;
    dv[1] = sys->statev[0].vec[1] - mouse_v.y;

    dxmag = mag(as_vec2(dx));

    float ks = 20;
    float kd = 1;

    f[0] = -(ks * (dxmag - 30) + kd * dv[0] * dx[0] / dxmag) * dx[0] / dxmag;
    f[1] = -(ks * (dxmag - 30) + kd * dv[1] * dx[1] / dxmag) * dx[1] / dxmag;

    sys->force[0].vec[0] += f[0];
    sys->force[0].vec[1] += f[1];
}


int
main(void)
{
    struct psysconfig config = {
        .boxh = sch,
        .boxw = scw,
        .cr = 0.9,
        .drag = DRAG,
        .gravity = GRAVITY,
        .invm = 1.0 / 10.0,
        .m = 10.0,
        .radius = 3.0,
    };
    struct psys sys;
    psys_init(&sys, config, 10);

    // da_append(&sys.forces, wall_collisions);
    da_append(&sys.forcecallbacks, mouse_coupling);
    da_append(&sys.forcecallbacks, spring);

    InitWindow(scw, sch, "physical-modelling-particles");
    SetTargetFPS(60);

    while (!WindowShouldClose()) {

        psys_step(&sys, delta_time);

        BeginDrawing();
        ClearBackground(BLACK);

        for (size_t i = 0; i < sys.count; i++) {
            for (size_t j = 0; j < sys.count; j++) {
                if (i == j) continue;
                DrawLineV(
                    *(Vector2 *)&sys.statex[i],
                    *(Vector2 *)&sys.statex[j],
                    (Color) {20, 20, 20, 255}
                );
            }
        }
        for (size_t i = 0; i < sys.count; i++) {
            // DrawPixelV(*(Vector2 *)sys.particles.items[i].x, RAYWHITE);
            DrawCircleV(*(Vector2 *)&sys.statex[i], sys.config.radius, RAYWHITE);
        }

        EndDrawing();
    }

    psys_delete(&sys);

    return 0;
}

