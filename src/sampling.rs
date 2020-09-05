extern crate nalgebra as na;

use std::cmp::max;

use na::geometry::Point2;
use na::Vector2;

use crate::pbrt::{Float, PI_OVER2, PI_OVER4};

use self::na::Vector3;

pub fn concentric_sample_disk(u: &Point2<Float>) -> Point2<Float> {
    let u_offset = (2.0 * u) - Vector2::new(1.0, 1.0);
    if u_offset.x == 0.0 && u_offset.y == 0.0 {
        return Point2::new(0.0, 0.0);
    }

    let (r, theta) = if u_offset.x.abs() > u_offset.y.abs() {
        (u_offset.x, PI_OVER4 * (u_offset.y / u_offset.x))
    } else {
        (u_offset.y, PI_OVER2 - PI_OVER4 * (u_offset.x / u_offset.y))
    };

    r * Point2::new(theta.cos(), theta.sin())
}

#[inline]
pub fn cosine_sample_hemisphere(u: &Point2<Float>) -> Vector3<Float> {
    let d = concentric_sample_disk(u);
    let sample = 1.0 - d.x * d.x - d.y * d.y;
    let z = sample.max(0.0).sqrt();
    Vector3::new(d.x, d.y, z)
}
