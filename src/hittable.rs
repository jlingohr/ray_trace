use super::point::Point;
use super::vector::Vec3;
use super::ray::Ray;
use super::material::Material;
use std::rc::Rc;


pub trait Hittable {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
}

#[derive(Clone)]
pub struct HitRecord {
    pub point: Point,
    pub normal: Vec3,
    pub t: f64,
    pub front_face: bool,
    pub material: Rc<dyn Material>,
}

impl HitRecord {
    pub fn new(point: Point, normal: Vec3, t: f64, front_face: bool, material: Rc<dyn Material>) -> HitRecord {
        HitRecord {
            point: point,
            normal: normal,
            t: t,
            front_face: front_face,
            material: material,
        }
    }
}

pub struct HittableList {
    pub objects: Vec<Rc<dyn Hittable>>,
}

impl HittableList {
    pub fn new() -> HittableList {
        HittableList { objects: Vec::new() }
    }

    pub fn clear(&mut self) {
        self.objects.clear();
    }

    pub fn add(&mut self, object: Rc<dyn Hittable>) {
        self.objects.push(object);
    }
}

impl Hittable for HittableList {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut hit_anything: Option<HitRecord> = None;
        let mut closest_so_far = t_max;

        for object in self.objects.iter() {
            if let Some(hit) = object.hit(ray, t_min, closest_so_far) {
                if hit.t < closest_so_far {
                    closest_so_far = hit.t;
                    hit_anything = Some(hit);
                }
            }
        }
        return hit_anything
    }
}