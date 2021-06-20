#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

const double PI = 3.1415926535897932385;
const double eps = 1e-8;

double degrees_to_radians(double degrees)
{
  return degrees * PI / 180.0;
}

double random_double()
{
  // real random in [0,1)
  return rand() / (RAND_MAX + 1.0);
}

double random_double_in_range(double min, double max) {
    // Returns a random real in [min,max).
    return min + (max-min)*random_double();
}

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

vec3 vec3_elementwise_multiply(vec3 a, vec3 b)
{
  vec3 ret;
  ret.x = a.x * b.x;
  ret.y = a.y * b.y;
  ret.z = a.z * b.z;
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

vec3 vec3_random_uniform(double min, double max)
{
  return vec3_new(
    random_double_in_range(min, max),
    random_double_in_range(min, max),
    random_double_in_range(min, max)
  );
}

vec3 vec3_random_in_unit_sphere()
{
  while(1)
  {
    vec3 p = vec3_random_uniform(-1.0, 1.0);
    if (vec3_length_squared(p) >= 1)
      continue;
    return p;
  }
}

vec3 vec3_random_unit_vector()
{
  return vec3_unit_vector(vec3_random_in_unit_sphere());
}

vec3 vec3_random_in_unit_hemisphere(vec3 normal)
{
  vec3 ret = vec3_random_in_unit_sphere();
  if (vec3_dot(ret, normal) > 0.0)
    return vec3_neg(ret);
  return ret;
}

int vec3_near_zero(vec3 v)
{
  return ((fabs(v.x) < eps) && (fabs(v.y) < eps) && (fabs(v.z) < eps));
}

vec3 reflect(vec3 vec, vec3 normal)
{
  return vec3_subtract(vec, vec3_scale(normal, 2*vec3_dot(vec, normal)));
}

vec3 refract(vec3 uv, vec3 n, double etai_over_etat)
{
  double cos_theta = fmin(vec3_dot(vec3_neg(uv), n), 1.0);
  vec3 r_out_perp = vec3_scale(vec3_add(uv, vec3_scale(n, cos_theta)), etai_over_etat);
  vec3 r_out_parl = vec3_scale(n, -sqrt(fabs(1.0-vec3_length_squared(r_out_perp))));
  return vec3_add(r_out_perp, r_out_parl);
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

ray ray_new(vec3 origin, vec3 direction)
{
  ray ret;
  ret.origin = origin;
  ret.direction = direction;
  return ret;
}

enum { MAT_LAMBERT, MAT_METAL, MAT_DIELECTRIC };

typedef struct material material;
struct material {
  int type;
  vec3 albedo;

  // metal
  double fuzz;

  // dielectric
  double ir;
};

typedef struct hit_record hit_record;
struct hit_record {
  vec3 p;
  vec3 normal;
  double t;
  int front_face; // are the normal and the ray acute ? 1:0
  material mat;
};

typedef struct sphere sphere;

struct sphere
{
  vec3 center;
  double radius;
  material mat;
};

int sphere_hit_set_face_normal(ray r, vec3 outward_normal, vec3* front_face)
{
  int is_front_face = vec3_dot(r.direction, outward_normal) < 0;
  *front_face = is_front_face ? outward_normal : vec3_neg(outward_normal);
  return is_front_face;
}

int hit_sphere(
  sphere s,
  double tmin,  double tmax,
  ray r,
  hit_record *hit
)
{
  vec3 oc = vec3_subtract(r.origin, s.center);
  double a = vec3_length_squared(r.direction);
  double half_b = vec3_dot(oc, r.direction);
  double c = vec3_length_squared(oc) - s.radius*s.radius;

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
  vec3 outward_normal = vec3_unit_vector(vec3_subtract(ray_point_at(r, root), s.center));
  hit->front_face = sphere_hit_set_face_normal(r,  outward_normal, &(hit->normal));
  hit->mat = s.mat;
  return 1;
}

double clamp(double x, double min, double max)
{
  if (x < min) return min;
  if (x > max) return max;
  return x;
}

void write_color_stdout(vec3 pixel_color, int samples_per_pixel)
{
  double r = pixel_color.x * (1.0 / samples_per_pixel);
  double g = pixel_color.y * (1.0 / samples_per_pixel);
  double b = pixel_color.z * (1.0 / samples_per_pixel);

  // apply gamma-2 and clamp between 0 and 256
  r = 256 * clamp(sqrt(r), 0.0, 0.99999);
  g = 256 * clamp(sqrt(g), 0.0, 0.99999);
  b = 256 * clamp(sqrt(b), 0.0, 0.99999);

  printf("%d %d %d\n", (int) r, (int) g, (int) b);
}

typedef struct world world;

struct world
{
  int n_objects;
  sphere* spheres;
};

double reflectance(double cosine, double ref_idx)
{
  double r0 = (1-ref_idx) / (1+ref_idx);
  r0 = r0*r0;
  return r0 + (1-r0)*pow((1 - cosine), 5);
}

int scatter(ray in_ray, vec3 at, vec3 normal, material mat, int front_face, ray* scattered_ray)
{
  int ret;
  switch (mat.type)
  {
    case MAT_LAMBERT:
      vec3 scatter_direction = vec3_add(normal, vec3_random_unit_vector());
      if (vec3_near_zero(scatter_direction))
      {
        scatter_direction = normal;
      }
      *scattered_ray = ray_new(at, scatter_direction);
      ret = 1;
      break;

    case MAT_METAL:
      vec3 reflected = reflect(vec3_unit_vector(in_ray.direction), normal);
      *scattered_ray = ray_new(at, vec3_add(reflected, vec3_scale(vec3_random_in_unit_sphere(), mat.fuzz)));
      ret = (vec3_dot(scattered_ray->direction, normal) > 0);
      break;

    case MAT_DIELECTRIC:
      double refraction_ratio = front_face ? (1.0/mat.ir) : mat.ir;
      vec3 unit_direction = vec3_unit_vector(in_ray.direction);
      double cos_theta = fmin(vec3_dot(vec3_neg(unit_direction), normal), 1.0);
      double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

      int cannot_refract = refraction_ratio * sin_theta > 1.0;

      vec3 new_ray;
      if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
        new_ray = reflect(unit_direction, normal);
      else
        new_ray = refract(unit_direction, normal, refraction_ratio);
      *scattered_ray = ray_new(at, new_ray);
      ret = 1;
  }
  return ret;
}

vec3 ray_color(ray r, world w, int depth)
{

  // max depth reached, no more rays
  if (depth <= 0) return vec3_new_zero();

  int ret = 0;
  hit_record rec;

  double min_t = 100000;
  hit_record min_rec;

  int hit = 0;

  for (int i = 0; i<w.n_objects; i++)
  {
    ret = hit_sphere(w.spheres[i], 0.001, DBL_MAX, r, &rec);
    if (ret)
    {
      hit = 1;
      if (rec.t < min_t)
      {
        min_t = rec.t;
        min_rec = rec;
      }
    }
  }

  if (hit == 1)
  {
    ray scattered_ray;
    int did_the_ray_scatter = scatter(r, min_rec.p, min_rec.normal, min_rec.mat, min_rec.front_face, &scattered_ray);
    if (did_the_ray_scatter)
    {
      return vec3_elementwise_multiply(
        ray_color(scattered_ray, w, depth-1),
      min_rec.mat.albedo);
    } else {
      return vec3_new(0.0, 0.0, 0.0);
    }
  }

  vec3 unit_direction = vec3_unit_vector(r.direction);
  double t = 0.5 * (unit_direction.y + 1.0);

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

typedef struct camera
{
  double aspect_ratio;
  int image_width, image_height;

  double viewport_height, viewport_width;
  double focal_length;

  vec3 origin;
  vec3 horizontal, vertical;
  vec3 lower_left_corner;
} camera;

camera camera_new_default()
{
  camera ret;
  ret.aspect_ratio = 16.0 / 9.0;
  ret.viewport_height = 2.0;
  ret.viewport_width = ret.aspect_ratio * ret.viewport_height;
  ret.focal_length = 1.0;

  ret.origin = vec3_new(0, 0, 0);
  ret.horizontal = vec3_new(ret.viewport_width, 0.0, 0.0);
  ret.vertical = vec3_new(0.0, ret.viewport_height, 0.0);

  vec3 origin_to_image_plane_center;
  origin_to_image_plane_center.x = 0;
  origin_to_image_plane_center.y = 0;
  origin_to_image_plane_center.z = ret.focal_length;

  // origin - 0.5*horizontal - 0.5*vertical - origin_to_image_plane_center
  ret.lower_left_corner = vec3_add(
    ret.origin, vec3_neg(
      vec3_add(
        vec3_add(
          vec3_scale(ret.horizontal, 0.5),
          vec3_scale(ret.vertical, 0.5)
        ),
        origin_to_image_plane_center
      )
    )
  );
  return ret;
}

ray camera_get_ray(camera cam, double u, double v)
{
  ray ret;
  ret.origin = cam.origin;

  // lower_left_corner + u*horizontal + v*vertical - origin
  ret.direction = vec3_add(
    vec3_add(
      vec3_add(
        vec3_scale(cam.horizontal, u),
        vec3_neg(cam.origin)
      ),
      vec3_scale(cam.vertical, v)
    ),
    cam.lower_left_corner
  );

  return ret;
}

int main()
{
  // image
  const double aspect_ratio = 16.0 / 9.0;
  const int image_width = 400;
  const int image_height = (int) (image_width / aspect_ratio);

  const int samples_per_pixel = 100;
  const int max_depth = 30;

  camera cam = camera_new_default();

  // world
  world w;
  w.n_objects = 4;


  // materials
  material mat_ground = {
    MAT_LAMBERT,
    vec3_new(0.8, 0.8, 0.0),
    0.0,
    0.0
  };

  material mat_center = {
    MAT_LAMBERT,
    vec3_new(0.1, 0.2, 0.5),
    0.0,
    1.5
  };

  material mat_left = {
    MAT_DIELECTRIC,
    vec3_new(1.0, 1.0, 1.0),
    0.3,
    1.2
  };

  material mat_right = {
    MAT_METAL,
    vec3_new(0.8, 0.6, 0.2),
    0.0,
    0.0
  };

  sphere ground =
  {
    vec3_new(0.0, -100.5, -1.0),
    100.0,
    mat_ground
  };

  sphere center =
  {
    vec3_new(0.0, 0.0, -1.0),
    0.5,
    mat_center
  };

  sphere left =
  {
    vec3_new(-1.0, 0.0, -1.0),
    0.5,
    mat_left
  };

  sphere right =
  {
    vec3_new(1.0, 0.0, -1.0),
    0.5,
    mat_right
  };

  w.spheres = (sphere[]) {ground, center, left, right};

  printf("P3\n%d %d\n255\n", image_width, image_height);

  for (int j = image_height-1; j>=0; --j)
  {
    fprintf(stderr, "Line %d of %d.\r", image_height-j, image_height);
    fflush(stderr);
    for (int i = 0; i<image_width; ++i)
    {
      vec3 px_color = vec3_new_zero();
      for (int s = 0; s < samples_per_pixel; s++)
      {
        double u = ((double) i + random_double()) / (image_width-1);
        double v = ((double) j + random_double()) / (image_height-1);
        ray r = camera_get_ray(cam, u, v);
        vec3 color = ray_color(r, w, max_depth);
        px_color = vec3_add(px_color, color);
      }
     write_color_stdout(px_color, samples_per_pixel);
     }
  }
  fprintf(stderr, "\nDone.\n");
}
