use super::aabb::AABB;
use super::material::Material;
use super::point::Point;
use super::ray::Ray;
use super::vector::Vec3;

pub trait Hittable {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB>;
}

#[derive(Clone)]
pub struct HitRecord<'a> {
    pub point: Point,
    pub normal: Vec3,
    pub t: f64,
    pub u: f64,
    pub v: f64,
    pub front_face: bool,
    pub material: &'a dyn Material,
}

impl<'a> HitRecord<'a> {
    pub fn new(
        point: Point,
        normal: Vec3,
        t: f64,
        u: f64,
        v: f64,
        front_face: bool,
        material: &'a dyn Material,
    ) -> HitRecord {
        HitRecord {
            point: point,
            normal: normal,
            t: t,
            u: u,
            v: v,
            front_face: front_face,
            material: material,
        }
    }
}

pub struct HittableList {
    pub objects: Vec<Box<dyn Hittable>>,
}

impl HittableList {
    pub fn new() -> HittableList {
        HittableList {
            objects: Vec::new(),
        }
    }

    pub fn clear(&mut self) {
        self.objects.clear();
    }

    pub fn add(&mut self, object: impl Hittable + 'static) {
        self.objects.push(Box::new(object));
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
        return hit_anything;
    }

    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        if self.objects.is_empty() {
            return None;
        }

        let mut output_box: Option<AABB> = None;

        // TODO this could probably be written in a fold
        for object in self.objects.iter() {
            match object.bounding_box(t0, t1) {
                None => return None,
                Some(bbox) => {
                    output_box = match output_box {
                        None => Some(bbox),
                        Some(temp_box) => Some(AABB::surrounding_box(&temp_box, &bbox)),
                    }
                }
            }
        }

        output_box
    }
}

pub struct FlipFace<H: Hittable> {
    pub hittable: H,
}

impl<H: Hittable> FlipFace<H> {
    pub fn new(hittable: H) -> FlipFace<H> {
        FlipFace { hittable }
    }
}

impl<H: Hittable> Hittable for FlipFace<H> {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        self.hittable.hit(ray, t_min, t_max).map(|mut hit| {
            hit.front_face = !hit.front_face;
            hit
        })
    }

    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        self.hittable.bounding_box(t0, t1)
    }
}

pub struct Translate<H: Hittable> {
    hittable: H,
    offset: Vec3,
}

impl<H: Hittable> Translate<H> {
    pub fn new(hittable: H, offset: Vec3) -> Translate<H> {
        Translate { hittable, offset }
    }
}

impl<H: Hittable> Hittable for Translate<H> {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let moved_ray = Ray::new(ray.origin - self.offset, ray.direction, ray.time);
        self.hittable.hit(&moved_ray, t_min, t_max).map(|mut hit| {
            hit.point += self.offset;
            hit.front_face = ray.direction.dot(&hit.normal) < 0.0;
            hit.normal = if hit.front_face {
                hit.normal
            } else {
                -hit.normal
            };
            hit
        })
    }

    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        self.hittable
            .bounding_box(t0, t1)
            .map(|bbox| AABB::new(bbox.min() + self.offset, bbox.max() + self.offset))
    }
}
