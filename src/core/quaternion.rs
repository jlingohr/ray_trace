use crate::core::pbrt::Float;

extern crate nalgebra as na;
use crate::core::transforms::Transform;
use na::Vector3;
use std::ops::{Add, Sub, Mul, Div, Neg, AddAssign, SubAssign, MulAssign, DivAssign};
use crate::core::math::clamp;

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Quaternion {
    pub v: Vector3<Float>,
    pub w: Float,
}

impl Quaternion {
    pub fn new(v: Vector3<Float>, w: Float) -> Quaternion {
        Quaternion {v, w,}
    }

    pub fn dot(&self, other: &Quaternion) -> Float {
        self.v.dot(&other.v) + (self.w * other.w)
    }

    pub fn normalize(&self) -> Quaternion {
        *self / self.dot(self).sqrt()
    }

    pub fn to_transform(&self) -> Transform {
        let xx = self.v.x * self.v.x;
        let yy = self.v.y * self.v.y;
        let zz = self.v.z * self.v.z;
        let xy = self.v.x * self.v.y;
        let xz = self.v.x * self.v.z;
        let yz = self.v.y * self.v.z;
        let wx = self.v.x * self.w;
        let wy = self.v.y * self.w;
        let wz = self.v.z * self.w;

        let mut mat = [[0.0; 4]; 4];
        mat[0][0] = 1.0 - (2.0 * (yy + zz));
        mat[0][1] = 2.0 * (xy + wz);
        mat[0][2] = 2.0 * (xz - wy);
        mat[1][0] = 2.0 * (xy - wz);
        mat[1][1] = 1.0 - (2.0 * (xx + zz));
        mat[1][2] = 2.0 * (yz + wx);
        mat[2][0] = 2.0 * (xz + wy);
        mat[2][1] = 2.0 * (yz - wx);
        mat[2][2] = 1.0 - (2.0 * (xx + yy));
        mat[3][3] = 1.0;

        Transform::new(mat).transpose()
    }

    pub fn from_transform(transform: &Transform) -> Quaternion {
        let m = &transform.matrix;
        let trace = m[(0, 0)] + m[(1, 1)] + m[(2, 2)];
        if trace > 0.0 {
            // compute w from matrix trace, then xyz
            let s = (trace + 1.0).sqrt();
            let w = s / 2.0;
            let s = 0.5 / s;
            let v = Vector3::new((m[(2, 1)] - m[(1, 2)]) * s, (m[(0, 2)] - m[(2, 0)]) * s, (m[(1, 0)] - m[(0, 1)]) * s);
            Quaternion {
                v,
                w
            }
        } else {
            // compute largest of x, y, z, then remaining components
            let nxt = [1, 2, 0];
            let mut q: [Float; 3] = [0.0; 3];
            let mut i = 0;
            if m[(1, 1)] > m[(0, 0)] {
                i = 1;
            }
            if m[(2, 2)] > m[(i, i)] {
                i = 2;
            }
            let j = nxt[i];
            let k = nxt[j];
            let mut s = m[(i, i)] - (m[(j, j)] + m[(k, k)] + 1.0).sqrt();
            q[i] = s * 0.5;
            if s != 0.0 {
                s = 0.5 / s;
            }
            let w = m[(k, j)] - m[(j, k)] * s;
            q[j] = (m[(j, i)] - m[(i, j)]) * s;
            q[k] = (m[(k, i)] + m[(i, k)]) * s;
            let v = Vector3::new(q[0], q[1], q[2]);
            Quaternion {
                v,
                w,
            }
        }
    }

    pub fn slerp(&self, other: &Quaternion, t: Float) -> Quaternion {
        let cos_theta = self.dot(other);
        if cos_theta > 0.9995 {
            (((1.0 - t) * self )+ (t * other)).normalize()
        } else {
            let theta = clamp(cos_theta, -1.0, 1.0).acos();
            let theta_p = theta * t;
            let qperp = (other - &(cos_theta * self)).normalize();
            theta_p.cos() * self + theta_p.sin() * qperp
        }
    }
}

impl AddAssign for Quaternion {
    fn add_assign(&mut self, other: Quaternion) {
        self.v += other.v;
        self.w += other.w;
    }
}

impl Add for Quaternion {
    type Output = Self;

    fn add(self, other: Quaternion) -> Quaternion {
        Quaternion {
            v: self.v + other.v,
            w: self.w + other.w,
        }
    }
}

impl SubAssign for Quaternion {
    fn sub_assign(&mut self, other: Quaternion) {
        self.v -= other.v;
        self.w -= other.w;
    }
}

impl Sub for Quaternion {
    type Output = Self;

    fn sub(self, other: Quaternion) -> Quaternion {
        Quaternion {
            v: self.v - other.v,
            w: self.w - other.w,
        }
    }
}

impl<'a, 'b> Sub<&'b Quaternion> for &'a Quaternion {
    type Output = Quaternion;

    fn sub(self, other: &'b Quaternion) -> Quaternion {
        Quaternion {
            v: &self.v - &other.v,
            w: &self.w - &other.w,
        }
    }
}

impl MulAssign<Float> for Quaternion {
    fn mul_assign(&mut self, other: Float) {
        self.v *= other;
        self.w += other;
    }
}

impl Mul<Float> for Quaternion {
    type Output = Self;

    fn mul(self, other: Float) -> Quaternion {
        Quaternion {
            v: self.v * other,
            w: self.w * other,
        }
    }
}

impl std::ops::Mul<Quaternion> for Float {
    type Output = Quaternion;

    fn mul(self, rhs: Quaternion) -> Quaternion {
        rhs * self
    }
}

impl<'a> std::ops::Mul<&'a Quaternion> for Float {
    type Output = Quaternion;

    fn mul(self, rhs: &Quaternion) -> Quaternion {
        Quaternion {
            v: &rhs.v * self,
            w: &rhs.w * self,
        }
    }
}

impl DivAssign<Float> for Quaternion {
    fn div_assign(&mut self, other: Float) {
        self.v /= other;
        self.w /= other;
    }
}

impl Div<Float> for Quaternion {
    type Output = Self;

    fn div(self, t: Float) -> Quaternion {
        Quaternion {
            v: self.v / t,
            w: self.w / t,
        }
    }
}

impl Neg for Quaternion {
    type Output = Self;

    fn neg(self) -> Quaternion {
        Quaternion {
            v: -self.v,
            w: -self.w,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    extern crate nalgebra as na;
    use na::Vector3;

    #[test]
    fn test_dot() {
        let q1 = Quaternion::new(Vector3::new(0.0, 1.0, 0.0), 1.0);
        let q2 = Quaternion::new(Vector3::new(1.0, 1.0, 1.0), 2.0);
        assert_eq!(3.0, q1.dot(&q2));
    }

    #[test]
    fn test_normalize() {
        let q1 = Quaternion::new(Vector3::new(1.0, 1.0, 1.0), 1.0);
        let expected = Quaternion::new(Vector3::new(0.5, 0.5, 0.5), 0.5);
        assert_eq!(expected, q1.normalize());
    }
}