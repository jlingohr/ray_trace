use super::aabb::AABB;
use super::hittable;
use super::material::Material;
use super::point::Point;
use super::ray::Ray;
use super::vector::Vec3;

#[derive(Clone)]
pub struct Sphere<M: Material> {
    pub center: Point,
    pub radius: f64,
    pub material: M,
}

impl<M: Material> Sphere<M> {
    pub fn new(center: Point, radius: f64, material: M) -> Sphere<M> {
        Sphere {
            center: center,
            radius: radius,
            material: material,
        }
    }
}

impl<M: Material> hittable::Hittable for Sphere<M> {
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
                let hit_record =
                    hittable::HitRecord::new(point, normal, temp, front_face, &self.material);
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
                let hit_record =
                    hittable::HitRecord::new(point, normal, temp, front_face, &self.material);
                return Some(hit_record);
            }
        }
        None
    }

    fn bounding_box(&self, _t0: f64, _t1: f64) -> Option<AABB> {
        let rad = Vec3::new(self.radius, self.radius, self.radius);
        let min = self.center - rad;
        let max = self.center + rad;
        let bbox = AABB::new(min, max);
        Some(bbox)
    }
}

pub struct MovingSphere<M: Material> {
    pub center0: Point,
    pub center1: Point,
    pub time0: f64,
    pub time1: f64,
    pub radius: f64,
    pub material: M,
}

impl<M: Material> MovingSphere<M> {
    pub fn new(
        center0: Point,
        center1: Point,
        time0: f64,
        time1: f64,
        radius: f64,
        material: M,
    ) -> MovingSphere<M> {
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

impl<M: Material> hittable::Hittable for MovingSphere<M> {
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
                let hit_record =
                    hittable::HitRecord::new(point, normal, temp, front_face, &self.material);
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
                let hit_record =
                    hittable::HitRecord::new(point, normal, temp, front_face, &self.material);
                return Some(hit_record);
            }
        }
        None
    }

    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        let rad = Vec3::new(self.radius, self.radius, self.radius);
        let box0 = AABB::new(self.center(t0) - rad, self.center(t0) + rad);
        let box1 = AABB::new(self.center(t1) - rad, self.center(t1) + rad);
        Some(AABB::surrounding_box(&box0, &box1))
    }
}
