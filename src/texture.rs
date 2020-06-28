use super::color::Color;
use super::perlin::PerlinNoise;
use super::point::Point;

pub trait Texture {
    fn value(&self, u: f64, v: f64, p: &Point) -> Color;
}

#[derive(Copy, Clone)]
pub struct SolidColor {
    color: Color,
}

impl SolidColor {
    pub fn new(r: f64, g: f64, b: f64) -> SolidColor {
        SolidColor {
            color: Color::new(r, g, b),
        }
    }
}

impl Texture for SolidColor {
    fn value(&self, _u: f64, _v: f64, _p: &Point) -> Color {
        self.color
    }
}

#[derive(Copy, Clone)]
pub struct Checkered<T: Texture> {
    odd: T,
    even: T,
}

impl<T: Texture> Checkered<T> {
    pub fn new(odd: T, even: T) -> Checkered<T> {
        Checkered { odd, even }
    }
}

impl<T: Texture> Texture for Checkered<T> {
    fn value(&self, u: f64, v: f64, p: &Point) -> Color {
        let sines = (10.0 * p.x).sin() * (10.0 * p.y).sin() * (10.0 * p.z).sin();
        if sines < 0.0 {
            self.odd.value(u, v, p)
        } else {
            self.even.value(u, v, p)
        }
    }
}

#[derive(Clone)]
pub struct NoiseTexture {
    noise: PerlinNoise,
    scale: f64,
}

impl NoiseTexture {
    pub fn new(scale: f64) -> NoiseTexture {
        let noise = PerlinNoise::new();
        NoiseTexture { noise, scale }
    }
}

impl Texture for NoiseTexture {
    fn value(&self, u: f64, v: f64, p: &Point) -> Color {
        let turb = self.noise.turb(p, 7);
        let amplitude = 1.0 + (self.scale * p.z + (10.0 * turb)).sin();
        Color::new(1.0, 1.0, 1.0) * 0.5 * amplitude
    }
}
