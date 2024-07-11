#include "raylib.h"
#include <math.h>
#include <raymath.h>
#include <stdio.h>

typedef struct {
  float r;
  Vector3 origin;
  Color color;
  int specular;
  float reflective;
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
Color trace_ray(Vector3 O, Vector3 D, float t_min, float t_max, size_t rec);
bool trace_shadow_ray(Vector3 P, Vector3 L, float t_min, float t_max);
void intersect_sphere(Vector3, Vector3, const sphere_t *, float *, float *);
float compute_lighting(Vector3 P, Vector3 N, Vector3 V, sphere_t *s);
Vector3 reflect_ray(Vector3 R, Vector3 N);

sphere_t spheres[4];
light_t lights[2];

int main(void) {
  InitWindow(cw, ch, "raylib [core] example - basic window");
  SetTargetFPS(60);
  Vector3 O = {.x = 0, .y = 0, .z = 0};

  spheres[0] = (sphere_t){.r = 1,
                          .origin = {.x = 0, .y = -1, .z = 3},
                          .color = RED,
                          .specular = 500,
                          .reflective = 0.2};
  spheres[1] = (sphere_t){.r = 1,
                          .origin = {.x = 2, .y = 0, .z = 4},
                          .color = BLUE,
                          .specular = 500,
                          .reflective = 0.3};
  spheres[2] = (sphere_t){.r = 1,
                          .origin = {.x = -2, .y = 0, .z = 4},
                          .color = GREEN,
                          .specular = 10,
                          .reflective = 0.4};
  spheres[3] = (sphere_t){.r = 5000,
                          .origin = {.x = 0, .y = -5001, .z = 0},
                          .color = YELLOW,
                          .specular = 1000,
                          .reflective = 0.5};
  lights[0] = (light_t){.type = POINT, .intensity = 0.6, .position = {2, 1, 0}};
  lights[1] =
      (light_t){.type = DIRECTIONAL, .intensity = 0.2, .position = {1, 4, 4}};

  while (!WindowShouldClose()) {
    BeginDrawing();
    ClearBackground(BLACK);
    for (int x = -cw / 2; x < cw / 2; x++) {
      for (int y = -ch / 2; y < ch / 2; y++) {
        Vector3 D = canvas_to_viewport(x, y);
        Color color = trace_ray(O, D, 1, INFINITY, 3);
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

Color trace_ray(Vector3 O, Vector3 D, float t_min, float t_max, size_t rec) {
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
    return BLACK;
  }
  Vector3 P = Vector3Add(O, Vector3Scale(D, closest_t));
  Vector3 ps = Vector3Subtract(P, closest_sphere.origin);
  Vector3 N = Vector3Normalize(ps);
  const Vector3 V = Vector3Scale(D, -1.0f);
  float i = compute_lighting(P, N, V, &closest_sphere);
  Color local_color = closest_sphere.color;
  local_color.r = Clamp(local_color.r * i, 0, 255);
  local_color.g = Clamp(local_color.g * i, 0, 255);
  local_color.b = Clamp(local_color.b * i, 0, 255);
  if (rec <= 0 || closest_sphere.reflective <= 0) {
    return local_color;
  }
  Vector3 R = reflect_ray(V, N);
  Color reflected_color = trace_ray(P, R, 0.1, INFINITY, rec - 1);
  reflected_color.r = reflected_color.r * closest_sphere.reflective +
                      local_color.r * (1 - closest_sphere.reflective);
  reflected_color.g = reflected_color.g * closest_sphere.reflective +
                      local_color.g * (1 - closest_sphere.reflective);
  reflected_color.b = reflected_color.b * closest_sphere.reflective +
                      local_color.b * (1 - closest_sphere.reflective);
  return reflected_color;
}
bool trace_shadow_ray(Vector3 P, Vector3 L, float t_min, float t_max) {

  for (size_t i = 0; i < 4; i++) {
    float t1, t2;
    intersect_sphere(P, L, &spheres[i], &t1, &t2);
    if ((t_min < t1 && t1 < t_max) || (t_min < t2 && t2 < t_max && t2)) {
      return true;
    }
  }
  return false;
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
    if (trace_shadow_ray(P, L, 0.001, INFINITY)) {
      continue;
    }
    float n_dot_l = Vector3DotProduct(N, L);
    if (n_dot_l > 0) {
      i += light->intensity * n_dot_l / (Vector3Length(N) * Vector3Length(L));
    }
    if (s->specular != -1) {
      Vector3 R = reflect_ray(L, N);
      float r_dot_v = Vector3DotProduct(R, V);
      if (r_dot_v > 0) {
        i += light->intensity *
             pow(r_dot_v / (Vector3Length(R) * Vector3Length(V)), s->specular);
      }
    }
  }
  return i;
}
Vector3 reflect_ray(Vector3 R, Vector3 N) {
  float n_dot_r = Vector3DotProduct(R, N);
  Vector3 ret = Vector3Scale(N, 2 * n_dot_r);
  ret = Vector3Subtract(ret, R);
  return ret;
}
