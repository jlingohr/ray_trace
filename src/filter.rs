use nalgebra::{Point2, Vector2};

use crate::pbrt::Float;

#[derive(Debug, Copy, Clone)]
pub enum Filter {
    Box(BoxFilter),
}

impl Filter {
    pub fn evaluate(&self, p: &Point2<Float>) -> Float {
        match self {
            Filter::Box(filter) => filter.evaluate(p),
        }
    }

    pub fn get_radius(&self) -> &Vector2<Float> {
        match self {
            Filter::Box(filter) => &filter.radius,
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub struct BoxFilter {
    pub radius: Vector2<Float>,
    pub inv_radius: Vector2<Float>,
}

impl BoxFilter {
    pub fn new(radius: Vector2<Float>) -> BoxFilter {
        let inv_radius = Vector2::new(1.0 / radius.x, 1.0 / radius.y);
        BoxFilter { radius, inv_radius }
    }

    fn evaluate(&self, p: &Point2<Float>) -> Float {
        1.0
    }
}
