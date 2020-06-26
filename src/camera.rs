use super::point::Point;
use super::vector::Vec3;
use super::ray::Ray;
use rand::prelude::ThreadRng;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Camera {
    origin: Point,
    horizontal: Vec3,
    vertical: Vec3,
    lower_left_corner: Point,
    u: Vec3,
    v: Vec3,
    w: Vec3,
    lens_radius: f64,
}

impl Camera {
    pub fn new(look_from: Point,
        look_at: Point,
        view_up: Vec3,
        fov: f64,
        aspect_ratio: f64,
        aperature: f64,
        focus_dist: f64) -> Camera {
        let theta = degrees_to_radians(fov);
        let h = (theta / 2.0).tan();
        let viewport_height = 2.0 * h;
        let viewport_width = aspect_ratio * viewport_height;

        let w = (look_from - look_at).unit();
        let u = view_up.cross(&w).unit();
        let v = w.cross(&u);

        let origin = look_from;
        let horizontal = focus_dist * viewport_width * u;
        let vertical = focus_dist * viewport_height * v;
        let lower_left_corner = origin - (horizontal / 2.0) - (vertical / 2.0) - (focus_dist * w);

        let lens_radius = aperature / 2.0;

        Camera {
            origin,
            horizontal,
            vertical,
            lower_left_corner,
            u,
            v,
            w,
            lens_radius,
        }
    }

    pub fn get_ray(&self, rng: &mut ThreadRng, s: f64, t: f64) -> Ray {
        let rd = self.lens_radius * Vec3::random_in_unit_disk(rng);
        let offset = self.u * rd.x + self.v * rd.y;
        Ray::new(self.origin + offset, self.lower_left_corner + (s * self.horizontal) + (t * self.vertical) - self.origin - offset)
    }
}

fn degrees_to_radians(degrees: f64) -> f64 {
    degrees * std::f64::consts::PI / 180.0
}