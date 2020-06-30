use super::aabb::AABB;
use super::hittable::{HitRecord, Hittable};
use super::material::Material;
use super::point::Point;
use super::ray::Ray;
use super::vector::Vec3;

pub enum Plane {
    XY,
    XZ,
    YZ,
}

impl Plane {
    pub fn get_axes(&self) -> (usize, usize, usize) {
        match *self {
            Plane::XY => (2, 0, 1),
            Plane::XZ => (1, 0, 2),
            Plane::YZ => (0, 1, 2),
        }
    }
}

pub struct AARectangle<M: Material> {
    mp: M,
    a0: f64,
    a1: f64,
    b0: f64,
    b1: f64,
    k: f64,
    plane: Plane,
}

impl<M: Material> AARectangle<M> {
    pub fn new(a0: f64, a1: f64, b0: f64, b1: f64, k: f64, mp: M, plane: Plane) -> AARectangle<M> {
        AARectangle {
            mp,
            a0,
            a1,
            b0,
            b1,
            k,
            plane,
        }
    }
}

impl<M: Material> Hittable for AARectangle<M> {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let (k_axis, a_axis, b_axis) = self.plane.get_axes();
        let t = (self.k - ray.origin[k_axis]) / ray.direction[k_axis];
        if t < t_min || t > t_max {
            return None;
        }

        let a = ray.origin[a_axis] + (t * ray.direction[a_axis]);
        let b = ray.origin[b_axis] + (t * ray.direction[b_axis]);

        if a < self.a0 || a > self.a1 || b < self.b0 || b > self.b1 {
            return None;
        }

        let u = (a - self.a0) / (self.a1 - self.a0);
        let v = (b - self.b0) / (self.b1 - self.b0);
        let mut outward_normal = Vec3::new(0.0, 0.0, 0.0);
        outward_normal[k_axis] = 1.0;
        let front_face = ray.direction.dot(&outward_normal) < 0.0;
        let normal = if front_face {
            outward_normal
        } else {
            -outward_normal
        };
        Some(HitRecord::new(
            ray.at(t),
            normal,
            t,
            u,
            v,
            front_face,
            &self.mp,
        ))
    }

    fn bounding_box(&self, _t0: f64, _t1: f64) -> Option<AABB> {
        let bbox = AABB::new(
            Point::new(self.a0, self.b0, self.k - 0.0001),
            Point::new(self.a1, self.b1, self.k + 0.0001),
        );
        Some(bbox)
    }
}
