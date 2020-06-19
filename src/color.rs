use std::ops::{Add, Mul, Sub, Div, AddAssign, MulAssign, SubAssign, DivAssign};
use super::vector::Vec3;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Color {
    pub r: f64,
    pub g: f64,
    pub b: f64,
}

impl Color {
    pub fn new(r: f64, g: f64, b:f64) -> Color {
        Color {r, g, b}
    }

    pub fn write_color(&self) {
        let ir: u32  = (255.999 * self.r) as u32;
        let ig: u32 = (255.999 * self.g) as u32;
        let ib: u32 = (255.999 * self.b) as u32;
        println!("{} {} {}", ir, ig, ib)
    }
}

impl AddAssign for Color {
    fn add_assign(&mut self, other: Color) {
        self.r += other.r;
        self.g += other.g;
        self.b += other.b;
    }
}

impl Add for Color {
    type Output = Self;

    fn add(self, other: Color) -> Color {
        Color {
            r: self.r + other.r,
            g: self.g + other.g,
            b: self.b + other.b
        }
    }
}

impl Add<Vec3> for Color {
    type Output = Self;

    fn add(self, other: Vec3) -> Color {
        Color {
            r: self.r + other.x,
            g: self.g + other.y,
            b: self.b + other.z
        }
    }
}

impl Add<Color> for Vec3 {
    type Output = Color;

    fn add(self, other: Color) -> Color {
        Color {
            r: self.x + other.r,
            g: self.y + other.g,
            b: self.z + other.b
        }
    }
}

impl SubAssign for Color {
    fn sub_assign(&mut self, other: Color) {
        self.r -= other.r;
        self.g -= other.g;
        self.b -= other.b;
    }
}

impl Sub for Color {
    type Output = Self;

    fn sub(self, other: Color) -> Color {
        Color {
            r: self.r - other.r,
            g: self.g - other.g,
            b: self.b - other.b
        }
    }
}

impl MulAssign for Color {
    fn mul_assign(&mut self, other: Color) {
        self.r *= other.r;
        self.g *= other.g;
        self.b *= other.b;
    }
}

impl Mul for Color {
    type Output = Self;

    fn mul(mut self, other: Color) -> Color {
        self *= other;
        self
    }
}

impl Mul<f64> for Color {
    type Output = Self;

    fn mul(self, t: f64) -> Color {
        Color {
            r: self.r * t,
            g: self.g * t,
            b: self.b * t,
        }
    }
}

impl std::ops::Mul<Color> for f64 {
    type Output = Color;

    fn mul(self, rhs: Color) -> Color {
        rhs * self
    }
}

impl DivAssign for Color {
    fn div_assign(&mut self, other: Color) {
        self.r /= other.r;
        self.g /= other.g;
        self.b /= other.b;
    }
}

impl Div for Color {
    type Output = Self;

    fn div(mut self, other: Color) -> Color {
        self /= other;
        self
    }
}

impl Div<f64> for Color {
    type Output = Self;

    fn div(self, t: f64) -> Color {
        Color {
            r: self.r / t,
            g: self.g / t,
            b: self.b / t,
        }
    }
}