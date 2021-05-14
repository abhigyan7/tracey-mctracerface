#include "vec.h"
#include "ray.h"
#include <stdio.h>

void write_color_stdout(vec3 pixel_color)
{
  int ir = (int) (255.999*pixel_color.x);
  int ig = (int) (255.999*pixel_color.y);
  int ib = (int) (255.999*pixel_color.z);

  printf("%d %d %d\n", ir, ig, ib);
}

double hit_sphere(vec3 center, double radius, ray r)
{
  vec3 oc = vec3_add(r.origin, vec3_neg(center));
  double a = vec3_length_squared(r.direction);
  double b = 2.0 * vec3_dot(oc, r.direction);
  double c = vec3_length_squared(oc) - radius*radius;

  double discriminant = b*b - 4*a*c;
  if (discriminant < 0)
  {
    return -1.0;
  } else {
    return (-b - sqrt(discriminant) ) / (2.0 * a);
  }
}

vec3 ray_color(ray r)
{

  vec3 c;
  c.x = 0; c.y = 0; c.z = -1;
  double radius = 0.5;

  double t = hit_sphere(c, radius, r);
  if (t > 0.0)
  {
    vec3 unit_vector = vec3_unit_vector(vec3_subtract(ray_point_at(r, t), c));
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

