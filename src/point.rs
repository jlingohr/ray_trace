use std::ops::{Add, Sub, AddAssign};
use super::vector::Vec3;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point {
    pub fn new(x: f64, y: f64, z:f64) -> Point {
        Point {x, y, z}
    }
}

impl AddAssign<Vec3> for Point {
    fn add_assign(&mut self, other: Vec3) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl Add<Vec3> for Point {
    type Output = Self;

    fn add(self, other: Vec3) -> Point {
        Point {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z
        }
    }
}

impl Add<Point> for Vec3 {
    type Output = Point;

    fn add(self, other: Point) -> Point {
        other + self
    }
}

impl Sub<Point> for Point {
    type Output = Vec3;

    fn sub(self, other: Point) -> Vec3 {
        Vec3::new(self.x-other.x, self.y-other.y, self.z-other.z)
    }
}

impl Sub<Vec3> for Point {
    type Output = Point;

    fn sub(self, other: Vec3) -> Point {
        Point::new(self.x-other.x, self.y-other.y, self.z-other.z)
    }
}

impl Sub<Point> for Vec3 {
    type Output = Point;

    fn sub(self, other: Point) -> Point {
        other - self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    mod tests_add {
        use super::*;

        #[test]
        fn add_point_vector() {
            let a = Point::new(0.0, 0.0, 0.0);
            let b = Vec3::ones();
            let c = Point::new(1.0, 1.0, 1.0);
            assert_eq!(a+b, c);
        }

        #[test]
        fn add_vector_point() {
            let a = Point::new(0.0, 0.0, 0.0);
            let b = Vec3::ones();
            let c = Point::new(1.0, 1.0, 1.0);
            assert_eq!(b+a, c);
        }
    }

    mod tests_sub {
        use super::*;

        #[test]
        fn sub_point_vector() {
            let a = Point::new(0.0, 0.0, 0.0);
            let b = Vec3::ones();
            let c = Point::new(-1.0, -1.0, -1.0);
            assert_eq!(a-b, c);
        }

        #[test]
        fn sub_vector_point() {
            let a = Point::new(0.0, 0.0, 0.0);
            let b = Vec3::ones();
            let c = Point::new(-1.0, -1.0, -1.0);
            assert_eq!(b-a, c);
        }

        #[test]
        fn sub_point_point() {
            let a = Point::new(0.0, 0.0, 0.0);
            let b = Point::new(1.0, 1.0, 1.0);
            let c = Vec3::new(-1.0, -1.0, -1.0);
            assert_eq!(a-b, c);
        }
    }
}
