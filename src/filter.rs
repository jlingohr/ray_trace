use crate::core::pbrt::Float;
use nalgebra::{Point2, Vector2};

pub trait Filter {
    fn evaluate(&self, p: &Point2<Float>) -> Float;
}

struct BoxFilter {
    pub radius: Vector2<Float>,
    pub inv_radius: Vector2<Float>,
}

impl BoxFilter {
    pub fn new(radius: Vector2<Float>) -> BoxFilter {
        let inv_radius = Vector2::new(1.0 / radius.x, 1.0 / radius.y);
        BoxFilter { radius, inv_radius }
    }
}

impl Filter for BoxFilter {
    fn evaluate(&self, p: &Point2<Float>) -> Float {
        1.0
    }
}
