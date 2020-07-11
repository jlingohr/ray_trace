extern crate nalgebra as na;
use nalgebra::geometry::{Point3, Point4};
use std::ops::Mul;

// TODO nalgebra uses column-major order so should optimize for this
// TODO handle error properly

pub struct Transform {
    matrix: na::Matrix4<f64>,
    matrix_inv: na::Matrix4<f64>,
}

impl Transform {
    pub fn new(mat: [[f64; 4]; 4]) -> Transform {
        let matrix = na::Matrix4::from(mat);
        if let Some(matrix_inv) = matrix.try_inverse() {
            Transform { matrix, matrix_inv }
        } else {
            panic!("Transformation matrix not invertible");
        }
    }

    pub fn inverse(t: &Transform) -> Transform {
        Transform {
            matrix: t.matrix_inv,
            matrix_inv: t.matrix,
        }
    }

    pub fn transpose(t: &Transform) -> Transform {
        Transform {
            matrix: t.matrix.transpose(),
            matrix_inv: t.matrix_inv.transpose(),
        }
    }

    pub fn translate(delta: &na::Vector3<f64>) -> Transform {
        let matrix = na::Matrix4::new(
            1.0, 0.0, 0.0, delta.x, 0.0, 1.0, 0.0, delta.y, 0.0, 0.0, 1.0, delta.z, 0.0, 0.0, 0.0,
            1.0,
        );
        let matrix_inv = na::Matrix4::new(
            1.0, 0.0, 0.0, -delta.x, 0.0, 1.0, 0.0, -delta.y, 0.0, 0.0, 1.0, -delta.z, 0.0, 0.0,
            0.0, 1.0,
        );
        Transform { matrix, matrix_inv }
    }

    pub fn scale(x: f64, y: f64, z: f64) -> Transform {
        let matrix = na::Matrix4::new(
            x, 0.0, 0.0, 0.0, 0.0, y, 0.0, 0.0, 0.0, 0.0, z, 0.0, 0.0, 0.0, 0.0, 1.0,
        );
        let matrix_inv = na::Matrix4::new(
            1.0 / x,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0 / y,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0 / z,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
        );
        Transform { matrix, matrix_inv }
    }

    pub fn rotate_x(theta: f64) -> Transform {
        let sin_theta = theta.sin();
        let cos_theta = theta.cos();
        let matrix = na::Matrix4::new(
            1.0, 0.0, 0.0, 0.0, 0.0, cos_theta, -sin_theta, 0.0, 0.0, sin_theta, cos_theta, 0.0,
            0.0, 0.0, 0.0, 1.0,
        );
        Transform {
            matrix: matrix,
            matrix_inv: matrix.transpose(),
        }
    }

    pub fn rotate_y(theta: f64) -> Transform {
        let sin_theta = theta.sin();
        let cos_theta = theta.cos();
        let matrix = na::Matrix4::new(
            cos_theta, 0.0, sin_theta, 0.0, 0.0, 1.0, 0.0, 0.0, -sin_theta, 0.0, cos_theta, 0.0,
            0.0, 0.0, 0.0, 1.0,
        );
        Transform {
            matrix: matrix,
            matrix_inv: matrix.transpose(),
        }
    }

    pub fn rotate_z(theta: f64) -> Transform {
        let sin_theta = theta.sin();
        let cos_theta = theta.cos();
        let matrix = na::Matrix4::new(
            cos_theta, -sin_theta, 0.0, 0.0, sin_theta, cos_theta, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0,
        );
        Transform {
            matrix: matrix,
            matrix_inv: matrix.transpose(),
        }
    }

    pub fn rotate(theta: f64, axis: &na::Vector3<f64>) -> Transform {
        let a = axis.normalize();
        let sin_theta = theta.sin();
        let cos_theta = theta.cos();
        let mut matrix = na::Matrix4::identity();
        // Rotation of first basis vector
        matrix[(0, 0)] = a.x * a.x + (1.0 - a.x * a.x) * cos_theta;
        matrix[(0, 1)] = a.x * a.x * (1.0 - cos_theta) - a.z * sin_theta;
        matrix[(0, 2)] = a.x * a.x * (1.0 - cos_theta) + a.y * sin_theta;
        matrix[(0, 3)] = 0.0;
        // Rotation of second basis vector
        matrix[(1, 0)] = a.x * a.y + (1.0 - cos_theta) + a.z * sin_theta;
        matrix[(1, 1)] = a.y * a.y + (1.0 - a.y * a.y) * cos_theta;
        matrix[(1, 2)] = a.y * a.z * (1.0 - cos_theta) - a.x * sin_theta;
        matrix[(1, 3)] = 0.0;
        // Rotation of third basis vector
        matrix[(2, 0)] = a.x * a.z * (1.0 - cos_theta) - a.y * sin_theta;
        matrix[(2, 1)] = a.y * a.z * (1.0 - cos_theta) + a.x * sin_theta;
        matrix[(2, 2)] = a.z * a.z + (1.0 - a.z * a.z) * cos_theta;
        matrix[(2, 3)] = 0.0;

        Transform {
            matrix: matrix,
            matrix_inv: matrix.transpose(),
        }
    }

    pub fn look_at(pos: &Point3<f64>, look: &Point3<f64>, up: &na::Vector3<f64>) -> Transform {
        let mut camera_to_world = na::Matrix4::identity();
        camera_to_world[(0, 3)] = pos.x;
        camera_to_world[(1, 3)] = pos.y;
        camera_to_world[(2, 3)] = pos.z;
        camera_to_world[(3, 3)] = 1.0;

        let dir = (look - pos).normalize();
        let right = (up.normalize().cross(&dir)).normalize();
        let new_up = dir.cross(&right);
        camera_to_world[(0, 0)] = right.x;
        camera_to_world[(1, 0)] = right.y;
        camera_to_world[(2, 0)] = right.z;
        camera_to_world[(3, 0)] = 0.0;
        camera_to_world[(0, 1)] = new_up.x;
        camera_to_world[(1, 1)] = new_up.y;
        camera_to_world[(2, 1)] = new_up.z;
        camera_to_world[(3, 1)] = 0.0;
        camera_to_world[(0, 2)] = dir.x;
        camera_to_world[(1, 2)] = dir.y;
        camera_to_world[(2, 2)] = dir.z;
        camera_to_world[(3, 2)] = 0.0;

        if let Some(matrix) = camera_to_world.try_inverse() {
            Transform {
                matrix: matrix,
                matrix_inv: camera_to_world,
            }
        } else {
            panic!("Camera to world matrix not invertible");
        }
    }

    pub fn transform_point(&self, p: &Point3<f64>) -> Point3<f64> {
        let p0 = Point4::new(p.x, p.y, p.z, 1.0);
        let p1 = self.matrix * p0;
        if p1.w == 1.0 {
            Point3::new(p1.x, p1.y, p1.z)
        } else {
            Point3::new(p1.x, p1.y, p1.z) / p1.w
        }
    }

    pub fn transform_vector(&self, v: &na::Vector3<f64>) -> na::Vector3<f64> {
        let v0 = na::Vector4::new(v.x, v.y, v.z, 0.0);
        let v1 = self.matrix * v0;
        na::Vector3::new(v1.x, v1.y, v1.z)
    }

    pub fn swap_handedness(&self) -> bool {
        self.matrix.determinant() < 0.0
    }
}

impl Mul for Transform {
    type Output = Self;

    fn mul(self, other: Transform) -> Transform {
        Transform {
            matrix: self.matrix * other.matrix,
            matrix_inv: other.matrix_inv * self.matrix_inv,
        }
    }
}
