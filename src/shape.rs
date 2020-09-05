use nalgebra::Point2;

use crate::geometry::bounds::Bounds3;
use crate::geometry::ray::Ray;
use crate::interaction::SurfaceInteraction;
use crate::material::Material;
use crate::pbrt::Float;
use crate::shapes::sphere::Sphere;

pub enum Shape {
    Sphere(Sphere),
}

impl Shape {
    pub fn object_bound(&self) -> Bounds3 {
        match self {
            Shape::Sphere(shape) => shape.object_bound(),
        }
    }

    pub fn world_bound(&self) -> Bounds3 {
        unimplemented!()
    }

    pub fn intersect(&self, ray: &Ray) -> Option<SurfaceInteraction> {
        match self {
            Shape::Sphere(shape) => shape.intersect(ray),
        }
    }
    pub fn intersect_p(&self, ray: &Ray) -> bool {
        let t_hit = ray.t_max;
        match self {
            Shape::Sphere(shape) => shape.intersect_p(ray),
        }
    }

    pub fn area(&self) -> Float {
        match self {
            Shape::Sphere(shape) => shape.area(),
        }
    }

    pub fn sample(&self, u: Point2<Float>) -> SurfaceInteraction {
        unimplemented!()
    }

    pub fn reverse_orientation(&self) -> bool {
        match self {
            Shape::Sphere(shape) => shape.reverse_orientation,
        }
    }

    pub fn get_material(&self) -> &Material {
        match self {
            Shape::Sphere(shape) => &shape.material,
        }
    }

    pub fn transform_swap_handedness(&self) -> bool {
        match self {
            Shape::Sphere(shape) => shape.transform_swap_handedness,
        }
    }
}
