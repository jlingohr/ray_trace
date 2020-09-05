use nalgebra::{Point3, Vector3};

use crate::pbrt::{next_float_down, next_float_up, Float};

pub mod bounds;
// mod lib;
pub mod ray;

pub fn offset_ray_origin(
    p: &Point3<Float>,
    p_error: &Vector3<Float>,
    n: &Vector3<Float>,
    w: &Vector3<Float>,
) -> Point3<Float> {
    let d = n.abs().dot(p_error);
    let mut offset = d * n;
    if w.dot(n) < 0.0 {
        offset = -offset;
    }
    let mut po = p + offset;
    for i in 0..3 {
        if offset[i] > 0.0 {
            po[i] = next_float_up(po[i as usize] as f32) as Float; //TODO should these be f64?
        } else if offset[i] < 0.0 {
            po[i] = next_float_down(po[i as usize] as f32) as Float;
        }
    }
    po
}

#[inline]
pub fn face_forward(a: Vector3<Float>, b: Vector3<Float>) -> Vector3<Float> {
    if a.dot(&b) < 0.0 {
        -a
    } else {
        a
    }
}
