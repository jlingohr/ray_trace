use super::hittable;
use super::material::Material;
use super::point::Point;
use super::ray::Ray;
use std::rc::Rc;

#[derive(Clone)]
pub struct Sphere {
    pub center: Point,
    pub radius: f64,
    pub material: Rc<dyn Material>,
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
                let normal = if front_face {
                    outward_normal
                } else {
                    -outward_normal
                };
                let hit_record = hittable::HitRecord::new(
                    point,
                    normal,
                    temp,
                    front_face,
                    Rc::clone(&self.material),
                );
                return Some(hit_record);
            }
            let temp = (-half_b + root) / a;
            if temp < t_max && temp > t_min {
                let point = ray.at(temp);
                let outward_normal = (point - self.center) / self.radius;
                let front_face = ray.direction.dot(&outward_normal) < 0.0;
                let normal = if front_face {
                    outward_normal
                } else {
                    -outward_normal
                };
                let hit_record = hittable::HitRecord::new(
                    point,
                    normal,
                    temp,
                    front_face,
                    Rc::clone(&self.material),
                );
                return Some(hit_record);
            }
        }
        None
    }
}

pub struct MovingSphere {
    pub center0: Point,
    pub center1: Point,
    pub time0: f64,
    pub time1: f64,
    pub radius: f64,
    pub material: Rc<dyn Material>,
}

impl MovingSphere {
    pub fn new(
        center0: Point,
        center1: Point,
        time0: f64,
        time1: f64,
        radius: f64,
        material: Rc<dyn Material>,
    ) -> MovingSphere {
        MovingSphere {
            center0,
            center1,
            time0,
            time1,
            radius,
            material,
        }
    }

    pub fn center(&self, time: f64) -> Point {
        self.center0
            + ((time - self.time0) / (self.time1 - self.time0)) * (self.center1 - self.center0)
    }
}

impl hittable::Hittable for MovingSphere {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<hittable::HitRecord> {
        let oc = ray.origin - self.center(ray.time);
        let a = ray.direction.dot(&ray.direction);
        let half_b = oc.dot(&ray.direction);
        let c = oc.dot(&oc) - (self.radius * self.radius);
        let discriminant = (half_b * half_b) - (a * c);

        if discriminant > 0.0 {
            let root = discriminant.sqrt();
            let temp = (-half_b - root) / a;
            if temp < t_max && temp > t_min {
                let point = ray.at(temp);
                let outward_normal = (point - self.center(ray.time)) / self.radius;
                let front_face = ray.direction.dot(&outward_normal) < 0.0;
                let normal = if front_face {
                    outward_normal
                } else {
                    -outward_normal
                };
                let hit_record = hittable::HitRecord::new(
                    point,
                    normal,
                    temp,
                    front_face,
                    Rc::clone(&self.material),
                );
                return Some(hit_record);
            }
            let temp = (-half_b + root) / a;
            if temp < t_max && temp > t_min {
                let point = ray.at(temp);
                let outward_normal = (point - self.center(ray.time)) / self.radius;
                let front_face = ray.direction.dot(&outward_normal) < 0.0;
                let normal = if front_face {
                    outward_normal
                } else {
                    -outward_normal
                };
                let hit_record = hittable::HitRecord::new(
                    point,
                    normal,
                    temp,
                    front_face,
                    Rc::clone(&self.material),
                );
                return Some(hit_record);
            }
        }
        None
    }
}
