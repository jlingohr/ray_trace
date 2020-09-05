extern crate nalgebra as na;

use na::geometry::Point3;
use na::Vector3;

use crate::pbrt::Float;

#[derive(Debug, Copy, Clone)]
pub struct Differential {
    pub rx_origin: Point3<Float>,
    pub ry_origin: Point3<Float>,
    pub rx_direction: Vector3<Float>,
    pub ry_direction: Vector3<Float>,
}

impl Differential {
    fn new(
        rx_origin: Point3<Float>,
        ry_origin: Point3<Float>,
        rx_direction: Vector3<Float>,
        ry_direction: Vector3<Float>,
    ) -> Differential {
        Differential {
            rx_origin,
            ry_origin,
            rx_direction,
            ry_direction,
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Ray {
    pub origin: Point3<Float>,
    pub direction: Vector3<Float>,
    pub time: Float,
    pub t_max: Float,
    pub differential: Option<Differential>,
}

impl Ray {
    pub fn new(origin: Point3<Float>, direction: Vector3<Float>, time: Float, t_max: Float) -> Ray {
        Ray {
            origin,
            direction,
            time,
            t_max,
            differential: None,
        }
    }

    pub fn new_with_differential(
        origin: Point3<Float>,
        direction: Vector3<Float>,
        time: Float,
        t_max: Float,
        rx_origin: Point3<Float>,
        ry_origin: Point3<Float>,
        rx_direction: Vector3<Float>,
        ry_direction: Vector3<Float>,
    ) -> Ray {
        let differential = Differential::new(rx_origin, ry_origin, rx_direction, ry_direction);
        Ray {
            origin,
            direction,
            time,
            t_max,
            differential: Some(differential),
        }
    }

    pub fn with_differential(
        self,
        rx_origin: Point3<Float>,
        ry_origin: Point3<Float>,
        rx_direction: Vector3<Float>,
        ry_direction: Vector3<Float>,
    ) -> Ray {
        let differential = Differential::new(rx_origin, ry_origin, rx_direction, ry_direction);
        Ray {
            differential: Some(differential),
            ..self
        }
    }

    pub fn at(&self, t: Float) -> Point3<Float> {
        &self.origin + (t * &self.direction)
    }

    pub fn scale_differential(self, s: Float) -> Ray {
        let scaled_differential = self.differential.map(|differential| {
            let rx_origin = &self.origin + (differential.rx_origin - &self.origin) * s;
            let ry_origin = &self.origin + (differential.ry_origin - &self.origin) * s;
            let rx_direction = &self.direction + (differential.rx_direction - &self.direction) * s;
            let ry_direction = &self.direction + (differential.ry_direction - &self.direction) * s;
            Differential::new(rx_origin, ry_origin, rx_direction, ry_direction)
        });

        Ray {
            differential: scaled_differential,
            ..self
        }
    }
}
