use crate::core::pbrt::{next_float_down, next_float_up, Float};
use std::ops::{Add, Div, Mul, Sub};

// Acts like a regular float but keeps track of an interval describing uncertainty of
// a value of interest, which arise due to errors in in floating-point arithmetic
#[derive(Debug, Copy, Clone)]
pub struct EFloat {
    v: f32,
    low: f32,
    high: f32,
}

impl EFloat {
    pub fn new(v: f32, err: f32) -> EFloat {
        if err == 0.0 {
            EFloat { v, low: v, high: v }
        } else {
            EFloat {
                v,
                low: next_float_down(v - err),
                high: next_float_up(v + err),
            }
        }
    }

    pub fn lower_bound(self) -> f32 {
        self.low
    }

    pub fn upper_bound(self) -> f32 {
        self.high
    }
}

impl PartialEq for EFloat {
    fn eq(&self, other: &Self) -> bool {
        self.v == other.v
    }
}

impl Add for EFloat {
    type Output = EFloat;

    fn add(self, rhs: Self) -> Self::Output {
        EFloat {
            v: self.v + rhs.v,
            low: next_float_down(self.lower_bound() + rhs.lower_bound()),
            high: next_float_up(self.upper_bound() + rhs.upper_bound()),
        }
    }
}

impl Sub for EFloat {
    type Output = EFloat;

    fn sub(self, rhs: Self) -> Self::Output {
        EFloat {
            v: self.v - rhs.v,
            low: next_float_down(self.lower_bound() - rhs.lower_bound()),
            high: next_float_up(self.upper_bound() - rhs.upper_bound()),
        }
    }
}

impl Mul for EFloat {
    type Output = EFloat;

    fn mul(self, rhs: Self) -> Self::Output {
        let prod: [f32; 4] = [
            self.lower_bound() * rhs.lower_bound(),
            self.upper_bound() * rhs.lower_bound(),
            self.lower_bound() * rhs.upper_bound(),
            self.upper_bound() * rhs.upper_bound(),
        ];

        EFloat {
            v: self.v * rhs.v,
            low: next_float_down(prod[0].min(prod[1]).min(prod[2].min(prod[3]))),
            high: next_float_up(prod[0].max(prod[1]).max(prod[2].max(prod[3]))),
        }
    }
}

impl Mul<f32> for EFloat {
    type Output = EFloat;
    fn mul(self, rhs: f32) -> EFloat {
        EFloat::new(rhs, 0.0) * self
    }
}

impl Div for EFloat {
    type Output = EFloat;

    fn div(self, rhs: Self) -> Self::Output {
        let div: [f32; 4] = [
            self.lower_bound() / rhs.lower_bound(),
            self.upper_bound() / rhs.lower_bound(),
            self.lower_bound() / rhs.upper_bound(),
            self.upper_bound() / rhs.upper_bound(),
        ];

        if rhs.low < 0.0 && rhs.high > 0.0 {
            EFloat {
                v: self.v / rhs.v,
                low: -std::f32::INFINITY,
                high: std::f32::INFINITY,
            }
        } else {
            EFloat {
                v: self.v / rhs.v,
                low: next_float_down(div[0].min(div[1]).min(div[2].min(div[3]))),
                high: next_float_up(div[0].max(div[1]).max(div[2].max(div[3]))),
            }
        }
    }
}
