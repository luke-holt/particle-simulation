#include <stdio.h>
#include <assert.h>

#include "raylib.h"

#define PASSERT UTIL_ASSERT
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

#define GRAVITY (-98.1)
#define DRAG (0.10)

const int fps = 60;
const float delta_time = 1.0 / fps;
const int scw = 800;
const int sch = 600;

void spring(pstate_t *sys) {
    for (size_t i = 0; i < sys->count; i++) {
        for (size_t j = 0; j < sys->count; j++) {
            if (i == j) continue;

            pvec_t dx = pvec_sub(sys->x[i], sys->x[j]);
            pvec_t dv = pvec_sub(sys->v[i], sys->v[j]);
            pvec_t dxnorm = pvec_normalize(dx);
            float dxmag = pvec_mag(dx);

            float ks = 1;
            float kd = 1;

            // -(ks * (dxmag - 30) + kd * dv * dx / dxmag) * dx / dxmag
            float a = ks * (dxmag - 30);
            float b = kd / dxmag;
            pvec_t f = {
                -(ks * (dxmag - 30) + kd * dv.vec[0] * dxnorm.vec[0]) * dxnorm.vec[0],
                -(ks * (dxmag - 30) + kd * dv.vec[1] * dxnorm.vec[1]) * dxnorm.vec[1],
            };
            pvec_accum(&sys->f[i], f);
            pvec_accum(&sys->f[j], pvec_scale(f, -1.0));
        }
    }
}

void mouse_coupling(pstate_t *sys) {
    if (!IsMouseButtonDown(MOUSE_BUTTON_LEFT))
        return;

    Vector2 mouse_x = GetMousePosition();
    Vector2 mouse_v = GetMouseDelta();

    pvec_t dx = pvec_sub(sys->x[0], *(pvec_t *)&mouse_x);
    pvec_t dv = pvec_sub(sys->v[0], *(pvec_t *)&mouse_v);
    pvec_t dxnorm = pvec_normalize(dx);
    float dxmag = pvec_mag(dx);

    float ks = 20;
    float kd = 1;

    // -(ks * (dxmag - 30) + kd * dv * dx / dxmag) * dx / dxmag
    float a = ks * (dxmag - 30);
    float b = kd / dxmag;
    pvec_t f = {
        -(ks * (dxmag - 30) + kd * dv.vec[0] * dxnorm.vec[0]) * dxnorm.vec[0],
        -(ks * (dxmag - 30) + kd * dv.vec[1] * dxnorm.vec[1]) * dxnorm.vec[1],
    };
    pvec_accum(&sys->f[0], f);
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
        .radius = 5,
    };
    pstate_t sys;
    particles_init(&sys, config, 10);

    particles_register_cb(&sys, mouse_coupling);
    // particles_register_cb(&sys, spring);

    InitWindow(scw, sch, "physical-modelling-particles");
    SetTargetFPS(60);

    while (!WindowShouldClose()) {

        particles_step(&sys, delta_time);

        BeginDrawing();
        ClearBackground(BLACK);

#if 0
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
#endif

        for (size_t i = 0; i < sys.count; i++) {
            // DrawPixelV(*(Vector2 *)sys.particles.items[i].x, RAYWHITE);
            DrawCircleV(*(Vector2 *)&sys.x[i], sys.config.radius, RAYWHITE);
        }

        EndDrawing();
    }

    particles_delete(&sys);

    return 0;
}

