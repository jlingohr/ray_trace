use super::aabb::AABB;
use super::aarect::{AARectangle, Plane};
use super::hittable::{FlipFace, HitRecord, Hittable, HittableList};
use super::material::Material;
use super::point::Point;
use super::ray::Ray;

pub struct Cube {
    box_min: Point,
    box_max: Point,
    sides: HittableList,
}

impl Cube {
    pub fn new<M: Material + Clone + 'static>(p0: Point, p1: Point, material: M) -> Cube {
        let box_min = p0;
        let box_max = p1;
        let mut sides = HittableList::new();
        sides.add(AARectangle::new(
            p0.x,
            p1.x,
            p0.y,
            p1.y,
            p1.z,
            material.clone(),
            Plane::XY,
        ));
        sides.add(FlipFace::new(AARectangle::new(
            p0.x,
            p1.x,
            p0.y,
            p1.y,
            p0.z,
            material.clone(),
            Plane::XY,
        )));
        sides.add(AARectangle::new(
            p0.x,
            p1.x,
            p0.z,
            p1.z,
            p1.y,
            material.clone(),
            Plane::XZ,
        ));
        sides.add(FlipFace::new(AARectangle::new(
            p0.x,
            p1.x,
            p0.z,
            p1.z,
            p0.y,
            material.clone(),
            Plane::XZ,
        )));
        sides.add(AARectangle::new(
            p0.y,
            p1.y,
            p0.z,
            p1.z,
            p1.x,
            material.clone(),
            Plane::YZ,
        ));
        sides.add(FlipFace::new(AARectangle::new(
            p0.y,
            p1.y,
            p0.z,
            p1.z,
            p0.x,
            material.clone(),
            Plane::YZ,
        )));

        Cube {
            box_min,
            box_max,
            sides,
        }
    }
}

impl Hittable for Cube {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        self.sides.hit(ray, t_min, t_max)
    }

    fn bounding_box(&self, _t0: f64, _t1: f64) -> Option<AABB> {
        Some(AABB::new(self.box_min, self.box_max))
    }
}
