use super::vector::Vec3;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Sub};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point {
    pub fn new(x: f64, y: f64, z: f64) -> Point {
        Point { x, y, z }
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
            z: self.z + other.z,
        }
    }
}

impl Add<Point> for Vec3 {
    type Output = Point;

    fn add(self, other: Point) -> Point {
        other + self
    }
}

impl Add<f64> for Point {
    type Output = Self;

    fn add(self, t: f64) -> Point {
        Point {
            x: self.x + t,
            y: self.y + t,
            z: self.z + t,
        }
    }
}

impl std::ops::Add<Point> for f64 {
    type Output = Point;

    fn add(self, rhs: Point) -> Point {
        rhs + self
    }
}

impl Sub<Point> for Point {
    type Output = Vec3;

    fn sub(self, other: Point) -> Vec3 {
        Vec3::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl Sub<Vec3> for Point {
    type Output = Point;

    fn sub(self, other: Vec3) -> Point {
        Point::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl Sub<Point> for Vec3 {
    type Output = Point;

    fn sub(self, other: Point) -> Point {
        other - self
    }
}

impl MulAssign<f64> for Point {
    fn mul_assign(&mut self, t: f64) {
        self.x *= t;
        self.y *= t;
        self.z *= t;
    }
}

impl Mul<f64> for Point {
    type Output = Self;

    fn mul(self, t: f64) -> Point {
        Point {
            x: self.x * t,
            y: self.y * t,
            z: self.z * t,
        }
    }
}

impl std::ops::Mul<Point> for f64 {
    type Output = Point;

    fn mul(self, rhs: Point) -> Point {
        rhs * self
    }
}

impl Index<usize> for Point {
    type Output = f64;

    fn index(&self, index: usize) -> &f64 {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Invalid index"),
        }
    }
}

impl IndexMut<usize> for Point {
    fn index_mut(&mut self, index: usize) -> &mut f64 {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("Invalid index"),
        }
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
            assert_eq!(a + b, c);
        }

        #[test]
        fn add_vector_point() {
            let a = Point::new(0.0, 0.0, 0.0);
            let b = Vec3::ones();
            let c = Point::new(1.0, 1.0, 1.0);
            assert_eq!(b + a, c);
        }
    }

    mod tests_sub {
        use super::*;

        #[test]
        fn sub_point_vector() {
            let a = Point::new(0.0, 0.0, 0.0);
            let b = Vec3::ones();
            let c = Point::new(-1.0, -1.0, -1.0);
            assert_eq!(a - b, c);
        }

        #[test]
        fn sub_vector_point() {
            let a = Point::new(0.0, 0.0, 0.0);
            let b = Vec3::ones();
            let c = Point::new(-1.0, -1.0, -1.0);
            assert_eq!(b - a, c);
        }

        #[test]
        fn sub_point_point() {
            let a = Point::new(0.0, 0.0, 0.0);
            let b = Point::new(1.0, 1.0, 1.0);
            let c = Vec3::new(-1.0, -1.0, -1.0);
            assert_eq!(a - b, c);
        }
    }
}
