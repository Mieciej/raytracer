#include "raylib.h"
#include <math.h>
#include <raymath.h>
#include <stdio.h>

typedef struct {
  float r;
  Vector3 origin;
  Color color;
  int specular;
} sphere_t;

const float ambient_light = 0.2;
typedef enum { POINT, DIRECTIONAL } light_type_t;
typedef struct {
  light_type_t type;
  float intensity;
  Vector3 position;
} light_t;

const int cw = 1000;
const int ch = 1000;

Vector3 canvas_to_viewport(int x, int y);
Color trace_ray(Vector3 O, Vector3 D, float t_min, float t_max);
void intersect_sphere(Vector3, Vector3, const sphere_t *, float *, float *);
float compute_lighting(Vector3 P, Vector3 N, Vector3 V, sphere_t *s);

sphere_t spheres[4];
light_t lights[2];

int main(void) {
  InitWindow(cw, ch, "raylib [core] example - basic window");
  SetTargetFPS(60);
  Vector3 O = {.x = 0, .y = 0, .z = 0};

  spheres[0] = (sphere_t){.r = 1,
                          .origin = {.x = 0, .y = -1, .z = 3},
                          .color = RED,
                          .specular = 100};
  spheres[1] = (sphere_t){.r = 1,
                          .origin = {.x = 2, .y = 0, .z = 4},
                          .color = BLUE,
                          .specular = 100};
  spheres[2] = (sphere_t){.r = 1,
                          .origin = {.x = -2, .y = 0, .z = 4},
                          .color = GREEN,
                          .specular = 500};
  spheres[3] = (sphere_t){.r = 5000,
                          .origin = {.x = 0, .y = -5001, .z = 0},
                          .color = YELLOW,
                          .specular = -1};
  lights[0] = (light_t){.type = POINT, .intensity = 0.6, .position = {2, 1, 0}};
  lights[1] =
      (light_t){.type = DIRECTIONAL, .intensity = 0.2, .position = {1, 4, 4}};

  while (!WindowShouldClose()) {
    BeginDrawing();
    ClearBackground(RAYWHITE);
    for (int x = -cw / 2; x < cw / 2; x++) {
      for (int y = -ch / 2; y < ch / 2; y++) {
        Vector3 D = canvas_to_viewport(x, y);
        Color color = trace_ray(O, D, 1, INFINITY);
        DrawPixel(x + cw / 2, -y + ch / 2, color);
      }
    }
    EndDrawing();
  }

  CloseWindow();

  return 0;
}
Vector3 canvas_to_viewport(int x, int y) {
  return (Vector3){.x = x * (1.0 / cw), .y = y * (1.0 / ch), .z = 1};
}

Color trace_ray(Vector3 O, Vector3 D, float t_min, float t_max) {
  float closest_t = INFINITY;
  sphere_t closest_sphere = {0};
  for (size_t i = 0; i < 4; i++) {
    float t1, t2;
    intersect_sphere(O, D, &spheres[i], &t1, &t2);
    if (t_min < t1 && t1 < t_max && t1 < closest_t) {
      closest_t = t1;
      closest_sphere = spheres[i];
    }
    if (t_min < t2 && t2 < t_max && t2 < closest_t) {
      closest_t = t2;
      closest_sphere = spheres[i];
    }
  }
  if (closest_sphere.r == 0) {
    return RAYWHITE;
  }
  Vector3 P = Vector3Add(O, Vector3Scale(D, closest_t));
  Vector3 ps = Vector3Subtract(P, closest_sphere.origin);
  Vector3 N = Vector3Normalize(ps);
  const Vector3 V = Vector3Scale(D, -1.0f);
  float i = compute_lighting(P, N, V, &closest_sphere);
  Color ret = closest_sphere.color;
  ret.r =Clamp(ret.r *  i,0,255);
  ret.g =Clamp(ret.g *  i,0,255);
  ret.b =Clamp(ret.b *  i,0,255);
  return ret;
}

void intersect_sphere(Vector3 O, Vector3 D, const sphere_t *sphere,
                      float *out_t1, float *out_t2) {
  float r = sphere->r;
  Vector3 CO = Vector3Subtract(O, sphere->origin);
  float a = Vector3DotProduct(D, D);
  float b = 2 * Vector3DotProduct(CO, D);
  float c = Vector3DotProduct(CO, CO) - r * r;
  float delta = b * b - 4 * a * c;
  if (delta < 0) {
    *out_t1 = INFINITY;
    *out_t2 = INFINITY;
    return;
  }
  *out_t1 = (-b + sqrtf(delta)) / (2 * a);
  *out_t2 = (-b - sqrtf(delta)) / (2 * a);
}

float compute_lighting(Vector3 P, Vector3 N, Vector3 V, sphere_t *s) {
  float i = 0.0f;
  i += ambient_light;
  Vector3 L;
  for (size_t j = 0; j < 2; j++) {
    light_t *light = &lights[j];
    if (light->type == POINT) {
      L = Vector3Subtract(light->position, P);
    } else {
      L = light->position;
    }
    float n_dot_l = Vector3DotProduct(N, L);
    if (n_dot_l > 0) {
      i += light->intensity * n_dot_l / (Vector3Length(N) * Vector3Length(L));
    }
    if (s->specular != -1) {
      float n_dot_l = Vector3DotProduct(N,L);
      Vector3 R = Vector3Scale(N,n_dot_l*2);
      R = Vector3Subtract(R,L);
      float r_dot_v = Vector3DotProduct(R, V);
      if (r_dot_v > 0) {
        i += light->intensity * pow(r_dot_v / (Vector3Length(R) * Vector3Length(V))
          , s->specular);
      }
    }
  }
  return i;
}
