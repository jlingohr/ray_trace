use super::vector::Vec3;
use super::point::Point;

#[derive(Copy, Clone)]
pub struct Ray {
    pub origin: Point,
    pub direction: Vec3,
}

impl Ray {
    pub fn new(origin: Point, direction: Vec3) -> Ray {
        Ray { origin, direction }
    }

    pub fn at(self, t: f64) -> Point {
        self.origin + (t*self.direction)
    }
}