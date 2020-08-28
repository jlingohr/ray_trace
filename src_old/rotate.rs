use super::aabb::AABB;
use super::hittable::{HitRecord, Hittable};
use super::point::Point;
use super::ray::Ray;

pub struct RotateY<H: Hittable> {
    pub hittable: H,
    pub sin_theta: f64,
    pub cos_theta: f64,
    pub bbox: Option<AABB>,
}

impl<H: Hittable> RotateY<H> {
    pub fn new(hittable: H, angle: f64) -> RotateY<H> {
        let radians = angle * (std::f64::consts::PI / 180.0);
        let sin_theta = radians.sin();
        let cos_theta = radians.cos();

        let bbox = hittable.bounding_box(0.0, 1.0).map(|bbox| {
            let mut min = Point::new(std::f64::INFINITY, std::f64::INFINITY, std::f64::INFINITY);
            let mut max = Point::new(
                -std::f64::INFINITY,
                -std::f64::INFINITY,
                -std::f64::INFINITY,
            );

            for i in 0..2 {
                let i = i as f64;
                for j in 0..2 {
                    let j = j as f64;
                    for k in 0..2 {
                        let k = k as f64;
                        let x = i * bbox.max().x + (1.0 - i) * bbox.min().x;
                        let y = j * bbox.max().y + (1.0 - j) * bbox.min().y;
                        let z = k * bbox.max().z + (1.0 - k) * bbox.min().z;

                        let x_new = cos_theta * x + sin_theta * z;
                        let z_new = -sin_theta * x + cos_theta * z;

                        if x_new < min.x {
                            min.x = x_new
                        };
                        if z_new < min.z {
                            min.z = z_new
                        };
                        if y < min.y {
                            min.y = y
                        };

                        if x_new > max.x {
                            max.x = x_new
                        };
                        if z_new > max.z {
                            max.z = z_new
                        };
                        if y > max.y {
                            max.y = y
                        };
                    }
                }
            }
            AABB::new(min, max)
        });
        RotateY {
            hittable,
            sin_theta,
            cos_theta,
            bbox,
        }
    }
}

impl<H: Hittable> Hittable for RotateY<H> {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut origin = ray.origin;
        let mut direction = ray.direction;

        origin.x = (self.cos_theta * ray.origin.x) - (self.sin_theta * ray.origin.z);
        origin.z = (self.sin_theta * ray.origin.x) + (self.cos_theta * ray.origin.z);

        direction.x = (self.cos_theta * ray.direction.x) - (self.sin_theta * ray.direction.z);
        direction.z = (self.sin_theta * ray.direction.x) + (self.cos_theta * ray.direction.z);

        let rotated_ray = Ray::new(origin, direction, ray.time);

        self.hittable
            .hit(&rotated_ray, t_min, t_max)
            .map(|mut hit| {
                let mut p = hit.point;
                let mut normal = hit.normal;

                p.x = (self.cos_theta * hit.point.x) + (self.sin_theta * hit.point.z);
                p.z = -(self.sin_theta * hit.point.x) + (self.cos_theta * hit.point.z);

                normal.x = (self.cos_theta * hit.normal.x) + (self.sin_theta * hit.normal.z);
                normal.z = -(self.sin_theta * hit.normal.x) + (self.cos_theta * hit.normal.z);

                hit.point = p;
                hit.front_face = rotated_ray.direction.dot(&hit.normal) < 0.0;
                hit.normal = if hit.front_face { normal } else { -normal };
                hit
            })
    }

    fn bounding_box(&self, _t0: f64, _t1: f64) -> Option<AABB> {
        self.bbox
    }
}
