use crate::core::spectrum::RGBSpectrum;

// Contains all global functions and constants
// Should be included in all other source files
pub type Float = f64;
pub type Spectrum = RGBSpectrum;

pub const PI: Float = 3.14159265358979323846;
pub const PI_OVER2: Float = 1.57079632679489661923;
pub const PI_OVER4: Float = 0.78539816339744830961;

pub const SHADOW_EPSILON: Float = 0.0001;
pub const ONE_MINUS_EPSILON: Float = 1.0 - Float::EPSILON;

pub const MACHINE_EPSILON: Float = std::f32::EPSILON * 0.5;

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

pub fn find_interval<P>(size: usize, pred: P) -> usize
where
    P: Fn(usize) -> bool,
{
    let mut first = 0;
    let mut len = size;
    while len > 0 {
        let half = len >> 1;
        let middle = first + half;
        if pred(middle) {
            first = middle + 1;
            len -= half + 1;
        } else {
            len = half;
        }
    }
    clamp(first - 1, 0, size - 2)
}

#[inline]
pub fn next_float_up(v: f32) -> f32 {
    if v.is_infinite() && v >= 0.0 {
        v
    } else {
        let new_v = if v == -0.0 { 0.0 } else { v };
        let mut ui = v.to_bits();
        if new_v >= 0.0 {
            ui += 1;
        } else {
            ui -= 1;
        }
        ui.to_f64()
    }
}

#[inline]
pub fn next_float_down(v: f32) -> f32 {
    if v.is_infinite() && v <= 0.0 {
        v
    } else {
        let new_v = if v == 0.0 { -0.0 } else { v };
        let mut ui = v.to_bits();
        if new_v > 0.0 {
            ui -= 1;
        } else {
            ui += 1;
        }
        ui.to_f64()
    }
}

// Find solution to  the quadratic equation at^2 + bt + c = 0
// Return Option<Float, Float> if solution is found, else None
#[inline]
pub fn quadratic(a: Float, b: Float, c: Float) -> Option<(Float, Float)> {
    let discrim = (b * b) - 4.0 * (a * c);
    if discrim < 0.0 {
        None
    } else {
        let root_discrim = discrim.sqrt();
        let q = if b < 0.0 {
            -0.5 * (b - root_discrim)
        } else {
            -0.5 * (b + root_discrim)
        };
        let mut t0 = q / a;
        let mut t1 = c / q;
        if t0 > t1 {
            std::mem::swap(&mut t0, &mut t1);
        }
        Some((t0, t1))
    }
}
