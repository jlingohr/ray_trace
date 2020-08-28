extern crate nalgebra as na;
use na::geometry::Point2;
use na::Vector2;
use crate::core::pbrt::{PI_OVER4, PI_OVER2, Float};

pub fn concentric_sample_disk(u: &Point2<Float>) -> Point2<Float> {
    let u_offset = (2.0 * u) - Vector2::new(1.0, 1.0);
    if u_offset.x == 0.0 && u_offset.y == 0.0 {
        return Point2::new(0.0, 0.0);
    }

    let (r, theta) = if u_offset.x.abs() > u_offset.y.abs() {
        (u_offset.x, PI_OVER4 * (u_offset.y / u_offset.x))
    } else {
       ( u_offset.y, PI_OVER2 - PI_OVER4 * (u_offset.x / u_offset.y))
    };

    r * Point2::new(theta.cos(), theta.sin())
}