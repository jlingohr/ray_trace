use super::point::Point;
use super::vector::Vec3;

#[derive(Copy, Clone)]
pub struct Ray {
    pub origin: Point,
    pub direction: Vec3,
    pub time: f64,
}

impl Ray {
    pub fn new(origin: Point, direction: Vec3, time: f64) -> Ray {
        Ray {
            origin,
            direction,
            time,
        }
    }

    pub fn at(self, t: f64) -> Point {
        self.origin + (t * self.direction)
    }
}
