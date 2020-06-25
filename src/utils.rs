use rand::prelude::ThreadRng;
use rand::Rng;

pub fn random_double(rng: &mut ThreadRng) -> f64 {
    rng.gen_range(0.0, 1.0)
}

pub fn random_in_range(rng: &mut ThreadRng, min: f64, max: f64) -> f64 {
    min + (max - min) * random_double(rng)
}

pub fn clamp(x: f64, min: f64, max: f64) -> f64 {
    if x < min {
        return min
    }
    if x > max {
        return max
    }
    return x
}
