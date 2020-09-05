extern crate nalgebra as na;

use std::ops::Mul;

use na::{Vector3, Vector4};
use nalgebra::geometry::Point3;

use crate::geometry::face_forward;
use crate::geometry::ray::Ray;
use crate::interaction::{CommonInteraction, Shading, SurfaceInteraction};
use crate::pbrt::{gamma, lerp, radians, Float};
use crate::quaternion::{slerp, Quaternion};

type Matrix4 = na::Matrix<Float, na::U4, na::U4, na::ArrayStorage<Float, na::U4, na::U4>>;

// TODO nalgebra uses column-major order so should optimize for this
// TODO handle error properly

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Transform {
    pub matrix: Matrix4,
    pub matrix_inv: Matrix4,
}

impl Default for Transform {
    fn default() -> Self {
        Transform {
            matrix: Matrix4::identity(),
            matrix_inv: Matrix4::identity(),
        }
    }
}

impl Transform {
    pub fn new(mat: [[Float; 4]; 4]) -> Transform {
        let matrix = Matrix4::new(
            mat[0][0], mat[0][1], mat[0][2], mat[0][3], mat[1][0], mat[1][1], mat[1][2], mat[1][3],
            mat[2][0], mat[2][1], mat[2][2], mat[2][3], mat[3][0], mat[3][1], mat[3][2], mat[3][3],
        );
        if let Some(matrix_inv) = matrix.try_inverse() {
            Transform { matrix, matrix_inv }
        } else {
            panic!("Transformation matrix not invertible");
        }
    }

    pub fn from_matrix(matrix: Matrix4) -> Transform {
        if let Some(matrix_inv) = matrix.try_inverse() {
            Transform { matrix, matrix_inv }
        } else {
            panic!("Transformation matrix not invertivle");
        }
    }

    pub fn inverse(&self) -> Transform {
        Transform {
            matrix: self.matrix_inv,
            matrix_inv: self.matrix,
        }
    }

    pub fn transpose(&self) -> Transform {
        Transform {
            matrix: self.matrix.transpose(),
            matrix_inv: self.matrix_inv.transpose(),
        }
    }

    // Create translation matrix from given vector
    // Only applies to points and should leave vectors unchanged
    pub fn translate(delta: &na::Vector3<Float>) -> Transform {
        let matrix = Matrix4::new(
            1.0, 0.0, 0.0, delta.x, 0.0, 1.0, 0.0, delta.y, 0.0, 0.0, 1.0, delta.z, 0.0, 0.0, 0.0,
            1.0,
        );
        let matrix_inv = Matrix4::new(
            1.0, 0.0, 0.0, -delta.x, 0.0, 1.0, 0.0, -delta.y, 0.0, 0.0, 1.0, -delta.z, 0.0, 0.0,
            0.0, 1.0,
        );
        Transform { matrix, matrix_inv }
    }

    pub fn scale(x: Float, y: Float, z: Float) -> Transform {
        let matrix = Matrix4::new(
            x, 0.0, 0.0, 0.0, 0.0, y, 0.0, 0.0, 0.0, 0.0, z, 0.0, 0.0, 0.0, 0.0, 1.0,
        );
        let matrix_inv = Matrix4::new(
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

    // Rotation along x-axis by theta degrees
    pub fn rotate_x(theta: Float) -> Transform {
        let sin_theta = radians(theta).sin();
        let cos_theta = radians(theta).cos();
        let matrix = Matrix4::new(
            1.0, 0.0, 0.0, 0.0, 0.0, cos_theta, -sin_theta, 0.0, 0.0, sin_theta, cos_theta, 0.0,
            0.0, 0.0, 0.0, 1.0,
        );
        Transform {
            matrix,
            matrix_inv: matrix.transpose(),
        }
    }

    pub fn rotate_y(theta: Float) -> Transform {
        let sin_theta = radians(theta).sin();
        let cos_theta = radians(theta).cos();
        let matrix = Matrix4::new(
            cos_theta, 0.0, sin_theta, 0.0, 0.0, 1.0, 0.0, 0.0, -sin_theta, 0.0, cos_theta, 0.0,
            0.0, 0.0, 0.0, 1.0,
        );
        Transform {
            matrix,
            matrix_inv: matrix.transpose(),
        }
    }

    pub fn rotate_z(theta: Float) -> Transform {
        let sin_theta = radians(theta).sin();
        let cos_theta = radians(theta).cos();
        let matrix = Matrix4::new(
            cos_theta, -sin_theta, 0.0, 0.0, sin_theta, cos_theta, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0,
        );
        Transform {
            matrix,
            matrix_inv: matrix.transpose(),
        }
    }

    pub fn rotate(theta: Float, axis: &na::Vector3<Float>) -> Transform {
        let a = axis.normalize();
        let sin_theta = radians(theta).sin();
        let cos_theta = radians(theta).cos();
        let mut matrix: Matrix4 = Matrix4::identity();
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
            matrix,
            matrix_inv: matrix.transpose(),
        }
    }

    pub fn orthographic(near: Float, far: Float) -> Transform {
        Transform::scale(1.0, 1.0, 1.0 / (far - near))
            * Transform::translate(&na::Vector3::new(0.0, 0.0, -near))
    }

    // Describes a transformation from a left-handed viewing coordinate system
    // Arguments:
    // * pos: Point3 - Camera position
    // * look: Point3 - Point being looked at
    // * up: Vector3 - Upwards direction
    pub fn look_at(
        pos: &Point3<Float>,
        look: &Point3<Float>,
        up: &na::Vector3<Float>,
    ) -> Transform {
        let mut camera_to_world: Matrix4 = Matrix4::identity();
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
                matrix,
                matrix_inv: camera_to_world,
            }
        } else {
            panic!("Camera to world matrix not invertible");
        }
    }

    // Applies a transformation to a 3d point by first converting the point to homogenous coordinates
    pub fn transform_point(&self, p: &Point3<Float>) -> Point3<Float> {
        let p0 = Vector4::new(p.x, p.y, p.z, 1.0);
        let p1 = &self.matrix * p0;
        if p1.w == 1.0 {
            Point3::new(p1.x, p1.y, p1.z)
        } else {
            Point3::new(p1.x, p1.y, p1.z) / p1.w
        }
    }

    // Applies a transformation to a 3d vector
    pub fn transform_vector(&self, v: &na::Vector3<Float>) -> Vector3<Float> {
        let v0 = na::Vector4::new(v.x, v.y, v.z, 0.0);
        let v1 = &self.matrix * v0;
        Vector3::new(v1.x, v1.y, v1.z)
    }

    pub fn transform_ray(&self, ray: &Ray) -> Ray {
        let (mut origin, o_error) = self.transform_point_with_error(&ray.origin);
        let direction = self.transform_vector(&ray.direction);
        let length = direction.dot(&direction);
        let mut t_max = ray.t_max;
        if length > 0.0 {
            let dt = direction.abs().dot(&o_error) / length;
            origin += &direction * dt;
            t_max -= dt;
        }

        Ray::new(origin, direction, ray.time, t_max)
    }

    pub fn transform_surface_interaction(&self, isect: &mut SurfaceInteraction) {
        let (point, p_error) = self.transform_point_with_error(&isect.interaction.point);
        isect.interaction.normal = self.transform_vector(&isect.interaction.normal).normalize();
        isect.interaction.w0 = self.transform_vector(&isect.interaction.w0);
        isect.interaction.time = isect.interaction.time;
        // isect.uv = isect.uv;
        isect.dpdu = self.transform_vector(&isect.dpdu);
        isect.dpdv = self.transform_vector(&isect.dpdv);
        isect.dndu = self.transform_vector(&isect.dndu);
        isect.dndv = self.transform_vector(&isect.dndv);
        isect.dudx = isect.dudx;
        isect.dvdx = isect.dvdx;
        isect.dudy = isect.dudy;
        isect.dvdy = isect.dvdy;
        isect.dpdx = self.transform_vector(&isect.dpdx);
        isect.dpdy = self.transform_vector(&isect.dpdy);
        isect.shading.normal = self.transform_vector(&isect.shading.normal).normalize();
        isect.shading.dpdu = self.transform_vector(&isect.shading.dpdu);
        isect.shading.dpdv = self.transform_vector(&isect.shading.dpdv);
        isect.shading.dndu = self.transform_vector(&isect.shading.dndu);
        isect.shading.dndv = self.transform_vector(&isect.shading.dndv);
        isect.shading.normal = face_forward(
            isect.shading.normal.clone_owned(),
            isect.interaction.normal.clone_owned(),
        );
        isect.dist = isect.dist;
    }

    pub fn swap_handedness(&self) -> bool {
        self.matrix.determinant() < 0.0
    }

    fn transform_point_with_error(&self, point: &Point3<Float>) -> (Point3<Float>, Vector3<Float>) {
        let p = Vector4::new(point.x, point.y, point.z, 1.0); //point.to_homogenous();
        let transformed_point = &self.matrix * p; //.transform_point(p);
        let abs_point = (&self.matrix * p).abs(); //.transform_point(p).abs(); // Not completely correct
        let p_error = Vector3::new(
            abs_point.x / abs_point.w,
            abs_point.y / abs_point.w,
            abs_point.z / abs_point.w,
        ) * gamma(3);
        if transformed_point.w == 1.0 {
            (
                Point3::new(
                    transformed_point.x,
                    transformed_point.y,
                    transformed_point.z,
                ),
                p_error,
            )
        } else {
            let inv = 1.0 / transformed_point.w;
            (
                Point3::new(
                    inv * transformed_point.x,
                    inv * transformed_point.y,
                    inv * transformed_point.z,
                ),
                p_error,
            )
        }
    }

    pub fn transform_vector_with_error(
        &self,
        v: &Vector3<Float>,
    ) -> (Vector3<Float>, Vector3<Float>) {
        let gamma = gamma(3);
        // TODO is there way to rust do Ax?
        let abs_error_x = gamma
            * ((self.matrix[(0, 0)] * v.x).abs()
                + (self.matrix[(0, 1)] * v.y).abs()
                + (self.matrix[(0, 2)] * v.z).abs());
        let abs_error_y = gamma
            * ((self.matrix[(1, 0)] * v.x).abs()
                + (self.matrix[(1, 1)] * v.y).abs()
                + (self.matrix[(1, 2)] * v.z).abs());
        let abs_error_z = gamma
            * ((self.matrix[(2, 0)] * v.x).abs()
                + (self.matrix[(2, 1)] * v.y).abs()
                + (self.matrix[(2, 2)] * v.z).abs());

        let transformed_vector = Vector3::new(
            self.matrix[(0, 0)] * v.x + self.matrix[(0, 1)] * v.y + self.matrix[(0, 2)] * v.z,
            self.matrix[(1, 0)] * v.x + self.matrix[(1, 1)] * v.y + self.matrix[(1, 2)] * v.z,
            self.matrix[(2, 0)] * v.x + self.matrix[(2, 1)] * v.y + self.matrix[(2, 2)] * v.z,
        );
        let v_error = Vector3::new(abs_error_x, abs_error_y, abs_error_z);

        (transformed_vector, v_error)
    }

    pub fn transform_ray_with_error(&self, ray: &Ray) -> (Ray, Vector3<Float>, Vector3<Float>) {
        let (mut o, o_error) = self.transform_point_with_error(&ray.origin);
        let (d, d_error) = self.transform_vector_with_error(&ray.direction);
        let length_squared = d.dot(&d);
        if length_squared > 0.0 {
            let dt = &d.abs().dot(&o_error) / length_squared;
            o += &d * dt;
        }
        let new_ray = Ray::new(o, d, ray.time, ray.t_max);

        (new_ray, o_error, d_error)
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

// Implement keyfram interpolation by decomposing given composite transformation matrices
// into scaling, rotation, and translation components
pub struct AnimatedTransform {
    start_transform: Transform,
    end_transform: Transform,
    start_time: Float,
    end_time: Float,
    is_animated: bool,
    translations: [Vector3<Float>; 2],
    rotations: [Quaternion; 2],
    scales: [Matrix4; 2],
    has_rotation: bool,
}

impl AnimatedTransform {
    pub fn new(
        start_transform: Transform,
        end_transform: Transform,
        start_time: Float,
        end_time: Float,
    ) -> AnimatedTransform {
        let is_animated = start_transform != end_transform;
        let (t0, r0, s0) = AnimatedTransform::decompose(&start_transform.matrix);
        let (t1, mut r1, s1) = AnimatedTransform::decompose(&end_transform.matrix);
        if r0.dot(&r1) < 0.0 {
            r1 = -r1;
        }
        let translations = [t0, t1];
        let rotations = [r0, r1];
        let scales = [s0, s1];
        let has_rotation = r0.dot(&r1) < 0.9995;
        // TODO compute terms of motion derivative function for bounding boxes

        AnimatedTransform {
            start_transform,
            end_transform,
            start_time,
            end_time,
            is_animated,
            translations,
            rotations,
            scales,
            has_rotation,
        }
    }

    fn decompose(matrix: &Matrix4) -> (Vector3<Float>, Quaternion, Matrix4) {
        // extract translation from transformation matrix
        let translation = Vector3::new(matrix[(0, 3)], matrix[(1, 3)], matrix[(2, 3)]);

        // compute new transformation M without translation
        let mut M = matrix.clone_owned();
        for i in 0..3 {
            M[(i, 3)] = 0.0;
            M[(3, i)] = 0.0;
        }
        M[(3, 3)] = 1.0;

        // extract rotation from transformation matrix
        let mut count = 0;
        let mut rotation: Matrix4 = M.clone_owned(); // should not need this?

        loop {
            // Compute next matrix Rnext
            // safe because matrix should already be validated on creation
            let rotation_it: Matrix4 = rotation.transpose().try_inverse().unwrap();
            let rotation_next: Matrix4 = 0.5 * (rotation + rotation_it);
            // Compute norm of different between R and RNext
            let norm: Float = (&rotation - &rotation_next).norm();
            rotation = rotation_next;
            count += 1;

            if count >= 100 || norm <= 0.0001 {
                break;
            }
        }
        let r_quat = Quaternion::from_transform(&Transform::from_matrix(rotation));
        // compute scale using rotation and original matrix
        let scale = rotation.try_inverse().unwrap() * M;

        return (translation, r_quat, scale);
    }

    pub fn transform_ray(&self, ray: &Ray) -> Ray {
        if !self.is_animated || ray.time <= self.start_time {
            self.start_transform.transform_ray(ray)
        } else if ray.time >= self.end_time {
            self.end_transform.transform_ray(ray)
        } else {
            let interpolation = self.interpolate(ray.time);
            interpolation.transform_ray(ray)
        }
    }

    pub fn transform_point(&self, time: Float, p: &Point3<Float>) -> Point3<Float> {
        if !self.is_animated || time < self.start_time {
            self.start_transform.transform_point(p)
        } else if time >= self.end_time {
            self.end_transform.transform_point(p)
        } else {
            let interpolation = self.interpolate(time);
            interpolation.transform_point(p)
        }
    }

    pub fn transform_vector(&self, time: Float, v: &Vector3<Float>) -> Vector3<Float> {
        if !self.is_animated || time < self.start_time {
            self.start_transform.transform_vector(v)
        } else if time >= self.end_time {
            self.end_transform.transform_vector(v)
        } else {
            let interpolation = self.interpolate(time);
            interpolation.transform_vector(v)
        }
    }

    fn interpolate(&self, time: Float) -> Transform {
        if !self.is_animated || time <= self.start_time {
            return self.start_transform;
        }
        if time >= self.end_time {
            return self.end_transform;
        }

        let dt = (time - self.start_time) / (self.end_time - self.start_time);

        // interpolate translation at dt
        let trans: Vector3<Float> = &self.translations[0] * (1.0 - dt) + dt * &self.translations[1];
        // interpolate rotation at dt
        let rotate = slerp(dt, &self.rotations[0], &self.rotations[1]);
        // interpolate scale at dt
        let mut scale = Matrix4::zeros();
        for i in 0..3 {
            for j in 0..3 {
                scale[(i, j)] = lerp(dt, self.scales[0][(i, j)], self.scales[1][(i, j)]);
            }
        }
        // compute interpolated matrix as product of interpolated components
        Transform::translate(&trans)
            * rotate.to_transform()
            * Transform {
                matrix: scale,
                matrix_inv: scale.try_inverse().unwrap(), //TODO better way
            }
    }
}

// Unit tests
#[cfg(test)]
mod test {
    use super::*;

    mod test_transform {
        use super::*;

        extern crate approx;

        #[test]
        fn test_transform_translate_point() {
            let transform = Transform::translate(&Vector3::new(1.0, 1.0, 1.0));
            let p_check = Point3::new(0.0, 0.0, 0.0);
            let p_expect = Point3::new(1.0, 1.0, 1.0);

            assert_eq!(p_expect, transform.transform_point(&p_check));
        }

        #[test]
        fn test_transform_translate_vector() {
            let transform = Transform::translate(&Vector3::new(1.0, 1.0, 1.0));
            let p_check = Vector3::new(0.0, 0.0, 0.0);
            let p_expect = Vector3::new(0.0, 0.0, 0.0);

            assert_eq!(p_expect, transform.transform_vector(&p_check));
        }

        #[test]
        fn test_transform_scale_point() {
            let transform = Transform::scale(0.5, 1.0, 2.0);
            let p_check = Point3::new(1.0, 1.0, 1.0);
            let p_expect = Point3::new(0.5, 1.0, 2.0);
            assert_eq!(p_expect, transform.transform_point(&p_check));
        }

        #[test]
        fn test_transform_scale_vector() {
            let transform = Transform::scale(0.5, 1.0, 2.0);
            let p_check = Vector3::new(1.0, 1.0, 1.0);
            let p_expect = Vector3::new(0.5, 1.0, 2.0);
            assert_eq!(p_expect, transform.transform_vector(&p_check));
        }

        #[test]
        fn test_transform_rotate_point_x() {
            let transform = Transform::rotate_x(90.0);
            let p_check = Point3::new(1.0, 1.0, 1.0);
            let p_expect = Point3::new(1.0, 0.0, -1.0);
            approx::abs_diff_eq!(p_expect, transform.transform_point(&p_check));
        }

        #[test]
        fn test_transform_rotate_point_y() {
            let transform = Transform::rotate_y(90.0);
            let p_check = Point3::new(1.0, 1.0, 1.0);
            let p_expect = Point3::new(-1.0, 1.0, 0.0);
            approx::abs_diff_eq!(p_expect, transform.transform_point(&p_check));
        }

        #[test]
        fn test_transform_rotate_point_z() {
            let transform = Transform::rotate_y(90.0);
            let p_check = Point3::new(1.0, 1.0, 1.0);
            let p_expect = Point3::new(-1.0, 0.0, 1.0);
            approx::abs_diff_eq!(p_expect, transform.transform_point(&p_check));
        }

        #[test]
        fn test_transform_rotate_point() {
            let transform = Transform::rotate(90.0, &Vector3::new(1.0, 1.0, 1.0));
            let p_check = Point3::new(1.0, 1.0, 1.0);
            let p_expect = Point3::new(0.0, 0.0, -1.0);
            approx::abs_diff_eq!(p_expect, transform.transform_point(&p_check));
        }
    }
}
