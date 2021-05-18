#include <stdio.h>
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

typedef struct ray
{
  vec3 origin;
  vec3 direction;
} ray;

vec3 ray_point_at(ray r, double t)
{
  return vec3_add(r.origin, vec3_scale(r.direction, t));
}

typedef struct hit_record hit_record;
struct hit_record {
  vec3 p;
  vec3 normal;
  double t;
};

int hit_sphere(
  vec3 center,  double radius,
  double tmin,  double tmax,
  ray r,
  hit_record *hit
)
{
  vec3 oc = vec3_subtract(r.origin, center);
  double a = vec3_length_squared(r.direction);
  double half_b = vec3_dot(oc, r.direction);
  double c = vec3_length_squared(oc) - radius*radius;

  double discriminant = half_b * half_b - a*c;
  if (discriminant < 0) return 0;
  double sqrtd = sqrt(discriminant);

  double root = (- half_b - sqrtd)/a;
  if (root < tmin || root > tmax)
  {
    root = (- half_b + sqrtd)/a;
    if (root < tmin || root > tmax)
      return 0;
  }
  hit->t = root;
  hit->p = ray_point_at(r, root);
  hit->normal = vec3_unit_vector(vec3_subtract(ray_point_at(r, root), center));
  return 1;
}
void write_color_stdout(vec3 pixel_color)
{
  int ir = (int) (255.999*pixel_color.x);
  int ig = (int) (255.999*pixel_color.y);
  int ib = (int) (255.999*pixel_color.z);

  printf("%d %d %d\n", ir, ig, ib);
}

vec3 ray_color(ray r)
{

  vec3 c;
  c.x = 0; c.y = 0; c.z = -1;
  double radius = 0.5;

  hit_record* hit = malloc(sizeof(hit_record));
  int t = hit_sphere(c, radius, 0, 1000, r, hit);
  if (t == 1)
  {
    vec3 unit_vector = hit->normal;
    vec3 color;
    color.x = unit_vector.x + 1.0;
    color.y = unit_vector.y + 1.0;
    color.z = unit_vector.z + 1.0;
    return vec3_scale(color, 0.5);
  }

  vec3 unit_direction = vec3_unit_vector(r.direction);
  t = 0.5 * (unit_direction.y + 1.0);

  vec3 color1;
  color1.x = 1.0; color1.y = 1.0; color1.z = 1.0;

  vec3 color2;
  color2.x = 0.5; color2.y = 0.7; color2.z = 1.0;

  t = 0.5 * (unit_direction.y + 1.0);
  return vec3_add(
    vec3_scale(color1, 1.0-t),
    vec3_scale(color2, t)
  );
}

int main()
{
  // image
  const double aspect_ratio = 16.0 / 9.0;
  const int image_width = 400;
  const int image_height = (int) (image_width / aspect_ratio);

  // camera
  
  double viewport_height = 2.0;
  double viewport_width = aspect_ratio * viewport_height;
  double focal_length = 1.0;

  vec3 origin;
  origin.x = 0;
  origin.y = 0;
  origin.z = 0;

  vec3 horizontal;
  horizontal.x = viewport_width; horizontal.y = 0; horizontal.z = 0;

  vec3 vertical;
  vertical.x = 0; vertical.y = viewport_height; vertical.z = 0;

  vec3 origin_to_image_plane_center;
  origin_to_image_plane_center.x = 0;
  origin_to_image_plane_center.y = 0;
  origin_to_image_plane_center.z = focal_length;

  // origin - 0.5*horizontal - 0.5*vertical - origin_to_image_plane_center
  vec3 lower_left_corner = vec3_add(
    origin, vec3_neg(
      vec3_add(
        vec3_add(
          vec3_scale(horizontal, 0.5),
          vec3_scale(vertical, 0.5)
        ), 
        origin_to_image_plane_center
      )
    )
  );

  // render
  printf("P3\n%d %d\n255\n", image_width, image_height);

  for (int j = image_height-1; j>=0; --j)
  {
    fprintf(stderr, "Line %d of %d.\r", image_height-j, image_height);
    fflush(stderr);
    for (int i = 0; i<image_width; ++i)
    {

      double u = ((double) i) / (image_width-1);
      double v = ((double) j) / (image_height-1);

      ray r;
      r.origin = origin;
      // lower_left_corner + u*horizontal + v*vertical - origin
      r.direction = vec3_add(
        vec3_add(
          vec3_add(
            vec3_scale(horizontal, u),
            vec3_neg(origin)
          ),
          vec3_scale(vertical, v)
        ),
        lower_left_corner
      );

     vec3 color = ray_color(r);
     write_color_stdout(color);
     }
  }
  fprintf(stderr, "\nDone.\n");
}

