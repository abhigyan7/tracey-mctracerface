#ifndef VEC3_H
#define VEC3_H

#include <math.h>
#include <stdlib.h>

typedef struct vec3
{
  double x;
  double y;
  double z;
} vec3;

vec3 vec3_new(double x, double y, double z)
{
  vec3 ret;
  ret.x = x;
  ret.y = y;
  ret.z = z;

  return ret;
}

vec3 vec3_new_zero()
{
  vec3 ret;
  ret.x = 0;
  ret.y = 0;
  ret.z = 0;
  return ret;
}

void vec3_del(vec3* vec)
{
  free(vec);
}

vec3 vec3_neg(vec3 vec)
{
  vec3 ret;
  ret.x = -vec.x;
  ret.y = -vec.y;
  ret.z = -vec.z;
  return ret;
}

double vec3_x(vec3 vec)
{
  return vec.x;
}

double vec3_y(vec3 vec)
{
  return vec.y;
}

double vec3_z(vec3 vec)
{
  return vec.z;
}	

vec3 vec3_add(vec3 a, vec3 b)
{
  vec3 ret;
  ret.x = a.x + b.x;
  ret.y = a.y + b.y;
  ret.z = a.z + b.z;
  return ret;
}

vec3 vec3_scale(vec3 a, double scalar)
{
  vec3 ret;
  ret.x = a.x * scalar;
  ret.y = a.y * scalar;
  ret.z = a.z * scalar; 
  return ret;
}

double vec3_length_squared(vec3 vec)
{
  return vec.x*vec.x + vec.y*vec.y + vec.z*vec.z;
}

double vec3_length(vec3 vec)
{
  return sqrt(vec3_length_squared(vec));
}

vec3 vec3_add_scalar(vec3 a, double scalar)
{
  vec3 ret;
  ret.x = a.x + scalar;
  ret.y = a.y + scalar;
  ret.z = a.z + scalar;
  return ret;
}

vec3 vec3_subtract(vec3 a, vec3 b)
{
  vec3 ret;
  ret.x = a.x - b.x;
  ret.y = a.y - b.y;
  ret.z = a.z - b.z;
  return ret;
}

vec3 vec3_cross(vec3 a, vec3 b)
{
  vec3 ret;
  ret.x = a.y * b.z - a.z * b.y;
  ret.y = a.z * b.x - a.x - b.z;
  ret.z = a.x * b.y - a.y * b.x;
  return ret;
}

double vec3_dot(vec3 a, vec3 b)
{
  return a.x * b.x 
       + a.y * b.y 
       + a.z * b.z;
}

vec3 vec3_unit_vector(vec3 v)
{
  vec3 ret;
  ret.x = v.x / vec3_length(v);
  ret.y = v.y / vec3_length(v);
  ret.z = v.z / vec3_length(v);
  return ret;
}

#endif

