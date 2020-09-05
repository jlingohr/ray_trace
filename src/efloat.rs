use std::ops::{Add, Div, Mul, Sub};

use crate::pbrt::{next_float_down, next_float_up, MACHINE_EPSILON};

// Acts like a regular float but keeps track of an interval describing uncertainty of
// a value of interest, which arise due to errors in in floating-point arithmetic
#[derive(Debug, Copy, Clone)]
pub struct EFloat {
    pub v: f32,
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

pub fn quadratic(a: EFloat, b: EFloat, c: EFloat) -> Option<(EFloat, EFloat)> {
    let discrim = (b.v as f64 * b.v as f64) - 4.0 * (a.v as f64 * c.v as f64);
    if discrim < 0.0 {
        None
    } else {
        let root_discrim = discrim.sqrt();
        let float_root_discrim = EFloat::new(
            root_discrim as f32,
            MACHINE_EPSILON as f32 * root_discrim as f32,
        );

        let q = if b.v < 0.0 {
            (b - float_root_discrim) * -0.5
        } else {
            (b + float_root_discrim) * -0.5
        };
        let mut t0 = q / a;
        let mut t1 = c / q;
        if t0.v > t1.v {
            std::mem::swap(&mut t0, &mut t1);
        }
        Some((t0, t1))
    }
}
