#ifndef RAY_H
#define RAY_H

#include "vec.h"

typedef struct ray
{
  vec3 origin;
  vec3 direction;
} ray;

vec3 ray_point_at(ray r, double t)
{
  return vec3_add(r.origin, vec3_scale(r.direction, t));
}

#endif

