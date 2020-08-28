use super::aabb::AABB;
use super::hittable::{HitRecord, Hittable};
use super::material::Isotropic;
use super::ray::Ray;
use super::texture::Texture;
use super::vector::Vec3;
use rand::Rng;

pub struct ConstantMedium<H: Hittable, T: Texture> {
    boundary: H,
    phase_function: Isotropic<T>,
    neg_inv_density: f64,
}

impl<H: Hittable, T: Texture> ConstantMedium<H, T> {
    pub fn new(boundary: H, d: f64, a: T) -> ConstantMedium<H, T> {
        let neg_inv_density = -1.0 / d;
        let phase_function = Isotropic::new(a);
        ConstantMedium {
            boundary,
            phase_function,
            neg_inv_density,
        }
    }
}

impl<H: Hittable, T: Texture> Hittable for ConstantMedium<H, T> {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        if let Some(mut rec1) = self
            .boundary
            .hit(ray, -std::f64::INFINITY, std::f64::INFINITY)
        {
            if let Some(mut rec2) = self.boundary.hit(ray, rec1.t + 0.0001, std::f64::INFINITY) {
                if rec1.t < t_min {
                    rec1.t = t_min;
                }
                if rec2.t > t_max {
                    rec2.t = t_max;
                }
                if rec1.t < rec2.t {
                    if rec1.t < 0.0 {
                        rec1.t = 0.0;
                    }
                    let ray_length = ray.direction.len();
                    let distance_inside_boundary = (rec2.t - rec1.t) * ray_length;
                    let mut rng = rand::thread_rng();
                    let hit_distance = self.neg_inv_density * rng.gen::<f64>().ln();

                    if hit_distance <= distance_inside_boundary {
                        let t = rec1.t + hit_distance / ray_length;
                        let p = ray.at(t);
                        let normal = Vec3::new(1.0, 0.0, 0.0);
                        let front_face = true;
                        return Some(HitRecord::new(
                            p,
                            normal,
                            t,
                            0.0,
                            0.0,
                            front_face,
                            &self.phase_function,
                        ));
                    }
                }
            }
        }
        None
    }

    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        self.boundary.bounding_box(t0, t1)
    }
}
