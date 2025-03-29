#include <stdio.h>
#include <assert.h>

#include "raylib.h"

#define PASSERT UTIL_ASSERT
#include "particles.h"

#define UTIL_IMPLEMENTATION
#include "util.h"

#include "vec.h"

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

void wall_collisions(void *ctx) {
    pstate_t *sys = (pstate_t *)ctx;
    for (size_t i = 0; i < sys->count; i++) {
        if ((sys->x[i].vec[0] - sys->config.radius) <= 0.0 && sys->v[i].vec[0] < 0.0 ||
            (sys->x[i].vec[0] + sys->config.radius) >= scw && sys->v[i].vec[0] > 0.0) {
            sys->v[i].vec[0] *= -1.0 * sys->config.cr;
        }
        if ((sys->x[i].vec[1] - sys->config.radius) <= 0.0 && sys->v[i].vec[1] < 0.0 ||
            (sys->x[i].vec[1] + sys->config.radius) >= sch && sys->v[i].vec[1] > 0.0) {
            sys->v[i].vec[1] *= -1.0 * sys->config.cr;
        }
    }
}

void spring(pstate_t *sys) {
    for (size_t i = 0; i < sys->count; i++) {
        for (size_t j = 0; j < sys->count; j++) {
            if (i == j) continue;

            float dx[2], dv[2], dxmag, f[2];
            dx[0] = sys->x[i].vec[0] - sys->x[j].vec[0];
            dx[1] = sys->x[i].vec[1] - sys->x[j].vec[1];
            dv[0] = sys->v[i].vec[0] - sys->v[j].vec[0];
            dv[1] = sys->v[i].vec[1] - sys->v[j].vec[1];

            dxmag = mag(as_vec2(dx));

            float ks = 1;
            float kd = 1;

            f[0] = -(ks * (dxmag - 30) + kd * dv[0] * dx[0] / dxmag) * dx[0] / dxmag;
            f[1] = -(ks * (dxmag - 30) + kd * dv[1] * dx[1] / dxmag) * dx[1] / dxmag;

            sys->f[i].vec[0] += f[0];
            sys->f[i].vec[1] += f[1];
            sys->f[j].vec[0] -= f[0];
            sys->f[j].vec[1] -= f[1];
        }
    }
}

void mouse_coupling(pstate_t *sys) {
    if (!IsMouseButtonDown(MOUSE_BUTTON_LEFT))
        return;

    Vector2 mouse_x = GetMousePosition();
    Vector2 mouse_v = GetMouseDelta();

    float dx[2], dv[2], dxmag, f[2];
    dx[0] = sys->x[0].vec[0] - mouse_x.x;
    dx[1] = sys->x[0].vec[1] - mouse_x.y;
    dv[0] = sys->v[0].vec[0] - mouse_v.x;
    dv[1] = sys->v[0].vec[1] - mouse_v.y;

    dxmag = mag(as_vec2(dx));

    float ks = 20;
    float kd = 1;

    f[0] = -(ks * (dxmag - 30) + kd * dv[0] * dx[0] / dxmag) * dx[0] / dxmag;
    f[1] = -(ks * (dxmag - 30) + kd * dv[1] * dx[1] / dxmag) * dx[1] / dxmag;

    sys->f[0].vec[0] += f[0];
    sys->f[0].vec[1] += f[1];
}


int
main(void)
{
    pconfig_t config = {
        .boxh = sch,
        .boxw = scw,
        .cr = 0.5,
        .drag = DRAG,
        .gravity = GRAVITY,
        .invm = 1.0 / 10.0,
        .m = 10.0,
        .radius = 3.0,
    };
    pstate_t sys;
    particles_init(&sys, config, 10);

    // particles_register_cb(&sys, wall_collisions);
    particles_register_cb(&sys, mouse_coupling);
    particles_register_cb(&sys, spring);

    InitWindow(scw, sch, "physical-modelling-particles");
    SetTargetFPS(60);

    while (!WindowShouldClose()) {

        particles_step(&sys, delta_time);

        BeginDrawing();
        ClearBackground(BLACK);

        for (size_t i = 0; i < sys.count; i++) {
            for (size_t j = 0; j < sys.count; j++) {
                if (i == j) continue;
                DrawLineV(
                    *(Vector2 *)&sys.x[i],
                    *(Vector2 *)&sys.x[j],
                    (Color) {20, 20, 20, 255}
                );
            }
        }
        for (size_t i = 0; i < sys.count; i++) {
            // DrawPixelV(*(Vector2 *)sys.particles.items[i].x, RAYWHITE);
            DrawCircleV(*(Vector2 *)&sys.x[i], sys.config.radius, RAYWHITE);
        }

        EndDrawing();
    }

    particles_delete(&sys);

    return 0;
}

