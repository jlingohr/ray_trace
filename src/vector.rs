use std::ops::{Add, Mul, Sub, Div, Neg, AddAssign, MulAssign, SubAssign, DivAssign};

use super::utils;
use rand::prelude::ThreadRng;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub fn new(x: f64, y: f64, z:f64) -> Vec3 {
        Vec3 {x, y, z}
    }

    pub fn ones() -> Vec3 {
        Vec3::new(1.0, 1.0, 1.0)
    }

    pub fn zeros() -> Vec3 {
        Vec3::new(0.0, 0.0, 0.0)
    }

    pub fn dot(&self, other: &Vec3) -> f64 {
        (self.x * other.x) + (self.y * other.y) + (self.z * other.z)
    }

    pub fn len(&self) -> f64 {
        self.dot(self).sqrt()
    }

    pub fn unit(&self) -> Vec3 {
        *self / self.len()
    }

    pub fn reflect(&self, other: Vec3) -> Vec3 {
        *self - 2.0 * self.dot(&other) * other
    }

    pub fn refract(&self, other: Vec3, etai_over_etat: f64) -> Vec3 {
        let cos_theta = -self.dot(&other);
        let r_out_parallel = etai_over_etat * (*self + (cos_theta*other));
        let r_out_perp = -(1.0 - (r_out_parallel.dot(&r_out_parallel))).sqrt() * other;
        r_out_parallel + r_out_perp
    }

    pub fn random(rng: &mut ThreadRng) -> Vec3 {
        Vec3::new(utils::random_double(rng), utils::random_double(rng), utils::random_double(rng))
    }

    pub fn random_unit_vector(rng: &mut ThreadRng) -> Vec3 {
        let a = utils::random_in_range(rng, 0.0, 2.0 * std::f64::consts::PI);
        let z = utils::random_in_range(rng, -1.0, 1.0);
        let r = (1.0 - z*z).sqrt();
        Vec3::new(r * a.cos(), r * a.sin(), z)
    }

    pub fn random_in_range(rng: &mut ThreadRng, min: f64, max: f64) -> Vec3 {
        Vec3::new(utils::random_in_range(rng, min, max), utils::random_in_range(rng, min, max), utils::random_in_range(rng, min, max))
    }

    pub fn random_in_unit_sphere(rng: &mut ThreadRng) -> Vec3 {
        let p = loop {
            let p = Vec3::random_in_range(rng, -1.0, 1.0);
            if p.dot(&p) >= 1.0 {
                break p;
            }
        };
        return p
    }

    pub fn random_in_hemisphere(rng: &mut ThreadRng, normal: &Vec3) -> Vec3 {
        let in_unit_sphere = Vec3::random_in_unit_sphere(rng);
        if in_unit_sphere.dot(normal) > 0.0 { // In the same hemisphere as the normal
            return in_unit_sphere
        } else {
            return -in_unit_sphere
        }
    }
}

impl AddAssign for Vec3 {
    fn add_assign(&mut self, other: Vec3) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl Add for Vec3 {
    type Output = Self;

    fn add(self, other: Vec3) -> Vec3 {
        Vec3 {
            x: self.x+other.x,
            y: self.y+other.y,
            z: self.z+other.z
        }
    }
}

impl Add<f64> for Vec3 {
    type Output = Self;

    fn add(self, t: f64) -> Vec3 {
        Vec3 {
            x: self.x + t,
            y: self.y + t,
            z: self.z + t,
        }
    }
}

impl std::ops::Add<Vec3> for f64 {
    type Output = Vec3;

    fn add(self, rhs: Vec3) -> Vec3 {
        rhs + self
    }
}

impl SubAssign for Vec3 {
    fn sub_assign(&mut self, other: Vec3) {
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
    }
}

impl Sub for Vec3 {
    type Output = Self;

    fn sub(self, other: Vec3) -> Vec3 {
        Vec3 {
            x: self.x-other.x,
            y: self.y-other.y,
            z: self.z-other.z
        }
    }
}

impl Sub<f64> for Vec3 {
    type Output = Self;

    fn sub(self, t: f64) -> Vec3 {
        Vec3 {
            x: self.x - t,
            y: self.y - t,
            z: self.z - t,
        }
    }
}

impl std::ops::Sub<Vec3> for f64 {
    type Output = Vec3;

    fn sub(self, rhs: Vec3) -> Vec3 {
        rhs - self
    }
}

impl MulAssign for Vec3 {
    fn mul_assign(&mut self, other: Vec3) {
        self.x *= other.x;
        self.y *= other.y;
        self.z *= other.z;
    }
}

impl Mul for Vec3 {
    type Output = Self;

    fn mul(self, other: Vec3) -> Vec3 {
        Vec3 {
            x: self.x*other.x,
            y: self.y*other.y,
            z: self.z*other.z
        }
    }
}

impl MulAssign<f64> for Vec3 {
    fn mul_assign(&mut self, t: f64) {
        self.x *= t;
        self.y *= t;
        self.z *= t;
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;

    fn mul(self, t: f64) -> Vec3 {
        Vec3 {
            x: self.x * t,
            y: self.y * t,
            z: self.z * t,
        }
    }
}

impl std::ops::Mul<Vec3> for f64 {
    type Output = Vec3;

    fn mul(self, rhs: Vec3) -> Vec3 {
        rhs * self
    }
}

impl DivAssign for Vec3 {
    fn div_assign(&mut self, other: Vec3) {
        self.x /= other.x;
        self.y /= other.y;
        self.z /= other.z;
    }
}

impl Div for Vec3 {
    type Output = Self;

    fn div(self, other: Vec3) -> Vec3 {
        Vec3 {
            x: self.x/other.x,
            y: self.y/other.y,
            z: self.z/other.z
        }
    }
}

impl Div<f64> for Vec3 {
    type Output = Self;

    fn div(self, t: f64) -> Vec3 {
        Vec3 {
            x: self.x / t,
            y: self.y / t,
            z: self.z / t,
        }
    }
}

impl Neg for Vec3 {
    type Output = Self;
    fn neg(self) -> Vec3 {
        Vec3 {
            x: -1.0*self.x,
            y: -1.0*self.y,
            z: -1.0*self.z
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    mod tests_dot {
        use super::*;

        #[test]
        fn dot() {
            let a = Vec3::ones();
            let b = Vec3::ones();
            let dot = a.dot(&b);
            assert_eq!(dot, 3.0);
        }
    }

    mod tests_len {
        use super::*;

        #[test]
        fn len() {
            let a = Vec3::new(1.0, 1.0, 1.0);
            let length = a.len();
            assert_eq!(length, (3.0 as f64).sqrt());
        }
    }

    mod tests_add {
        use super::*;

        #[test]
        fn vector_add() {
            let a = Vec3::ones();
            let b = Vec3::ones();
            let c = Vec3::new(2.0, 2.0, 2.0);
            assert_eq!(a+b, c);
        }

        #[test]
        fn broadcast_add() {
            let a = Vec3::ones();
            let b: f64 = 1.0;
            let c = Vec3::new(2.0, 2.0, 2.0);
            assert_eq!(a+b, c);
        }
    }

    mod tests_sub {
        use super::*;

        #[test]
        fn vector_sub() {
            let a = Vec3::ones();
            let b = Vec3::ones();
            let c = Vec3::new(0.0, 0.0, 0.0);
            assert_eq!(a-b, c);
        }

        #[test]
        fn broadcast_sub() {
            let a = Vec3::ones();
            let b: f64 = 1.0;
            let c = Vec3::new(0.0, 0.0, 0.0);
            assert_eq!(a-b, c);
        }
    }

    mod tests_mul {
        use super::*;

        #[test]
        fn scalar_mul() {
            let t = 2.0;
            let a = Vec3::new(0.0, 1.0, 2.0);
            let b = Vec3::new(0.0, 2.0, 4.0);
            assert_eq!((t as f64)*a, b);
        }

        #[test]
        fn elem_mul() {
            let a = Vec3::new(0.0, 1.0, 2.0);
            let b = Vec3::new(1.0, 1.0, 1.0);
            let c = Vec3::new(0.0, 1.0, 2.0);
            assert_eq!(a*b, c);
        }
    }

    mod tests_div {
        use super::*;

        #[test]
        fn scalar_div() {
            let t = 2.0;
            let a = Vec3::new(1.0, 2.0, 4.0);
            let b = Vec3::new(0.5, 1.0, 2.0);
            assert_eq!(a/(t as f64), b);
        }

        #[test]
        fn elem_div() {
            let a = Vec3::new(0.0, 1.0, 2.0);
            let b = Vec3::new(1.0, 2.0, 2.0);
            let c = Vec3::new(0.0, 0.5, 1.0);
            assert_eq!(a/b, c);
        }
    }

    mod tests_neg {
        use super::*;

        #[test]
        fn neg_pos() {
            let a = Vec3::new(0.0, 1.0, 2.0);
            let b = Vec3::new(0.0, -1.0, -2.0);
            assert_eq!(-a, b);
        }

        #[test]
        fn neg_neg() {
            let a = Vec3::new(0.0, -1.0, -2.0);
            let b = Vec3::new(0.0, 1.0, 2.0);
            assert_eq!(-a, b);
        }
    }
}