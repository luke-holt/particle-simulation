#include <stdio.h>
#include <assert.h>

#include "raylib.h"

#define PASSERT UTIL_ASSERT
#include "particles.h"

#define UTIL_IMPLEMENTATION
#include "util.h"

#define GRAVITY 0.0 // (-981.)
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
                -(ks * (dxmag - 30) + kd * dv.x * dxnorm.x) * dxnorm.x,
                -(ks * (dxmag - 30) + kd * dv.y * dxnorm.y) * dxnorm.y,
            };
            pvec_accum(&sys->f[i], f);
            pvec_accum(&sys->f[j], pvec_scale(f, -1.0));
        }
    }
}

void gravitate_to_center(pstate_t *sys) {
    for (size_t i = 0; i < sys->count; i++) {
        pvec_t dx = pvec_sub(sys->x[i], (pvec_t){scw/2, sch/2});
        pvec_t dv = pvec_sub(sys->v[i], (pvec_t){0, 0});
        pvec_t dxnorm = pvec_normalize(dx);
        float dxmag = pvec_mag(dx);

        float ks = 10;
        float kd = 1;

        // -(ks * (dxmag - 30) + kd * dv * dx / dxmag) * dx / dxmag
        float a = ks * (dxmag - 30);
        float b = kd / dxmag;
        pvec_t f = {
            -(ks * (dxmag - 30) + kd * dv.x * dxnorm.x) * dxnorm.x,
            -(ks * (dxmag - 30) + kd * dv.y * dxnorm.y) * dxnorm.y,
        };
        pvec_accum(&sys->f[i], f);
    }
}


int
main(void)
{
    pconfig_t config = {
        .box = (pvec_t){scw, sch},
        .cr = 0.5,
        .drag = DRAG,
        .gravity = (pvec_t){0.0, -GRAVITY},
        .invm = 1.0 / 10.0,
        .m = 10.0,
        .radius = 5,
    };
    pstate_t sys;
    particles_init(&sys, config, 100);

    particles_register_cb(&sys, gravitate_to_center);
    // particles_register_cb(&sys, spring);

    InitWindow(scw, sch, "physical-modelling-particles");
    SetTargetFPS(60);

    double new_particle_timer = 0.0;

    while (!WindowShouldClose()) {

        double t = GetTime();
        if (new_particle_timer < t) {
            new_particle_timer = t + 0.25;
            pvec_t x = { config.radius * 2, config.radius * 2 };
            pvec_t v = { 0.0, 0.0 };
            particles_add(&sys, x, v);
        }

        particles_step(&sys, GetFrameTime());

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

