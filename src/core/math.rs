use crate::core::pbrt::{Float, PI};

#[inline]
pub fn radians(deg: Float) -> Float {
    (PI / 180.0) * deg
}

#[inline]
pub fn degrees(rad: Float) -> Float {
    (180.0 / PI) * rad
}

#[inline]
pub fn lerp(t: Float, v1: Float, v2: Float) -> Float {
    (1.0 - t) * v1 + (t * v2)
}

#[inline]
pub fn clamp<T: PartialOrd>(val: T, low: T, high: T) -> T {
    if val < low {
        low
    } else if val > high {
        high
    } else {
        val
    }
}