use super::point::Point;
use super::ray::Ray;

#[derive(Clone, Copy)]
pub struct AABB {
    min: Point,
    max: Point,
}

impl AABB {
    pub fn new(a: Point, b: Point) -> AABB {
        AABB { min: a, max: b }
    }

    pub fn min(&self) -> Point {
        self.min
    }

    pub fn max(&self) -> Point {
        self.max
    }

    pub fn hit(&self, ray: &Ray, tmin: f64, tmax: f64) -> bool {
        for a in 0..3 {
            let inv_d = 1.0 / ray.direction[a];
            let mut t0 = (self.min[a] - ray.origin[a]) * inv_d;
            let mut t1 = (self.max[a] - ray.origin[a]) * inv_d;
            if inv_d < 0.0 {
                std::mem::swap(&mut t0, &mut t1);
            }
            let min = if t0 > tmin { t0 } else { tmin };
            let max = if t1 < tmax { t1 } else { tmax };
            if max <= min {
                return false;
            }
        }
        return true;
    }

    pub fn surrounding_box(box0: &AABB, box1: &AABB) -> AABB {
        let small = Point::new(
            f64::min(box0.min.x, box1.min.x),
            f64::min(box0.min.y, box1.min.y),
            f64::min(box0.min.z, box1.min.z),
        );
        let big = Point::new(
            f64::max(box0.max.x, box1.max.x),
            f64::max(box0.max.y, box1.max.y),
            f64::max(box0.max.z, box1.max.z),
        );
        AABB::new(small, big)
    }
}
