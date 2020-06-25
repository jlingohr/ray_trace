use super::hittable;
use super::point::Point;
use super::ray::Ray;
use super::material::Material;
use std::rc::Rc;

#[derive(Clone)]
pub struct Sphere {
    pub center: Point,
    pub radius: f64,
    pub material:  Rc<dyn Material>,
}

impl Sphere {
    pub fn new(center: Point, radius: f64, material: Rc<dyn Material>) -> Sphere {
        Sphere {
            center: center,
            radius: radius,
            material: material,
        }
    }
}

impl hittable::Hittable for Sphere {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<hittable::HitRecord> {
        let oc = ray.origin - self.center;
        let a = ray.direction.dot(&ray.direction);
        let half_b = oc.dot(&ray.direction);
        let c = oc.dot(&oc) - (self.radius * self.radius);
        let discriminant = (half_b * half_b) - (a * c);

        if discriminant > 0.0 {
            let root = discriminant.sqrt();
            let temp = (-half_b - root) / a;
            if temp < t_max && temp > t_min {
                let point = ray.at(temp);
                let outward_normal = (point - self.center) / self.radius;
                let front_face = ray.direction.dot(&outward_normal) < 0.0;
                let normal = if front_face { outward_normal } else { -outward_normal };
                let hit_record = hittable::HitRecord::new(point, normal, temp, front_face, Rc::clone(&self.material));
                return Some(hit_record)
            }
            let temp = (-half_b + root) / a;
            if temp < t_max && temp > t_min {
                let point = ray.at(temp);
                let outward_normal = (point - self.center) / self.radius;
                let front_face = ray.direction.dot(&outward_normal) < 0.0;
                let normal = if front_face { outward_normal } else { -outward_normal };
                let hit_record = hittable::HitRecord::new(point, normal, temp, front_face, Rc::clone(&self.material));
                return Some(hit_record)
            }
        }
        None
    }
}