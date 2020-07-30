extern crate nalgebra as na;
use crate::core::transforms::Transform;
use na::Vector3;
use std::ops::{Add, Sub, Mul, Div, AddAssign, SubAssign, MulAssign, DivAssign};

pub struct Quaternion {
    pub v: Vector3<f64>,
    pub w: f64,
}

impl Quaternion {
    pub fn new(v: Vector3<f64>, w: f64) -> Quaternion {
        Quaternion {v, w,}
    }

    pub fn dot(&self, other: &Quaternion) -> f64 {
        self.v.dot(&other.v) + (self.w * other.w)
    }

    pub fn normalize(&self) -> Quaternion {
        self / self.dot(self).sqrt()
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
        mat[(0, 0)] = 1.0 - (2.0 * (yy + zz));
        mat[(0, 1)] = 2.0 * (xy + wz);
        mat[(0, 2)] = 2.0 * (xz - wy);
        mat[(1, 0)] = 2.0 * (xy - wz);
        mat[(1, 1)] = 1.0 - (2.0 * (xx + zz));
        mat[(1, 2)] = 2.0 * (yz + wx);
        mat[(2, 0)] = 2.0 * (xz + wy);
        mat[(2, 1)] = 2.0 * (yz - wx);
        mat[(2, 2)] - 1.0 - (2.0 * (xx + yy));
        mat[(3, 3)] = 1.0;

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
            let q = [f64; 3];
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

    pub fn slerp(&self, other: &Quaternion, t: f64) -> Quaternion {
        let cos_theta = self.dot(other);
        if cos_theta > 0.9995 {
            ((1.0 - t) * self + t * other).normalize()
        } else {
            let theta = cos_theta.clamp(-1.0, 1.0).acos();
            let theta_p = theta * t;
            let qperp = (other - q1 * cos_theta).normalize();
            self * theta_p.cos() + qperp * theta_p.sin()
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

impl MulAssign<f64> for Quaternion {
    fn mul_assign(&mut self, other: f64) {
        self.v *= other;
        self.w += other;
    }
}

impl Mul<f64> for Quaternion {
    type Output = Self;

    fn mul(self, other: other) -> Quaternion {
        Quaternion {
            v: self.v * other,
            w: self.w * other,
        }
    }
}

impl DivAssign<f64> for Quaternion {
    fn div_assign(&mut self, other: f64) {
        self.v /= other;
        self.w /= other;
    }
}

impl Div<f64> for Quaternion {
    type Output = Self;

    fn div(self, other: other) -> Quaternion {
        Quaternion {
            v: self.v / other,
            w: self.w / other,
        }
    }
}