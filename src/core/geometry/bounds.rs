use crate::core::pbrt::Float;

extern crate nalgebra as na;
use na::geometry::{Point2, Point3};
use na::{Vector2, Vector3};

use std::ops::Index;

// TODO should be able to type parametrerize bounds2<T> but couldn't
// get to work with nalgebra
#[derive(Debug, Copy, Clone)]
pub struct Bounds2i {
    pub p_min: Point2<i32>,
    pub p_max: Point2<i32>,
}

impl Bounds2i {
    pub fn new(p0: Point2<i32>, p1: Point2<i32>) -> Bounds2i {
        let p_min = Point2::new(p0.x.min(p1.x), p0.y.min(p1.y));
        let p_max = Point2::new(p0.x.max(p1.x), p0.y.max(p1.y));
        Bounds2i { p_min, p_max }
    }

    pub fn corner(&self, corner: usize) -> Point2<i32> {
        let x = if corner & 1 == 0 {
            self.p_min.x
        } else {
            self.p_max.x
        };
        let y = if corner & 2 == 0 {
            self.p_min.y
        } else {
            self.p_max.y
        };

        Point2::new(x, y)
    }

    pub fn union(&self, p: &Point2<i32>) -> Bounds2i {
        let p0 = Point2::new(self.p_min.x.min(p.x), self.p_min.y.min(p.y));
        let p1 = Point2::new(self.p_max.x.max(p.x), self.p_max.y.max(p.y));
        Bounds2i {
            p_min: p0,
            p_max: p1,
        }
    }

    pub fn union_bounds(&self, b: &Bounds2i) -> Bounds2i {
        let p0 = Point2::new(self.p_min.x.min(b.p_min.x), self.p_min.y.min(b.p_min.y));
        let p1 = Point2::new(self.p_max.x.max(b.p_max.x), self.p_max.y.max(b.p_max.y));
        Bounds2i {
            p_min: p0,
            p_max: p1,
        }
    }

    pub fn intersect(&self, b: &Bounds2i) -> Bounds2i {
        let p0 = Point2::new(self.p_min.x.max(b.p_min.x), self.p_min.y.max(b.p_min.y));
        let p1 = Point2::new(self.p_max.x.min(b.p_max.x), self.p_max.y.min(b.p_max.y));
        Bounds2i {
            p_min: p0,
            p_max: p1,
        }
    }

    pub fn overlaps(&self, b: &Bounds2i) -> bool {
        (self.p_max >= b.p_min) && (self.p_min <= b.p_max)
    }

    pub fn inside(&self, p: &Point2<i32>) -> bool {
        (p >= &self.p_min) && (p <= &self.p_max)
    }

    pub fn inside_exclusive(&self, p: &Point2<i32>) -> bool {
        (p >= &self.p_min) && (p < &self.p_max)
    }

    pub fn expand(&self, delta: Float) -> Bounds2i {
        Bounds2i {
            p_min: &self.p_min - Vector2::new(delta, delta),
            p_max: &self.p_max + Vector2::new(delta, delta),
        }
    }

    pub fn diagonal(&self) -> Vector2<i32> {
        &self.p_max - &self.p_min
    }

    pub fn surface_area(&self) -> i32 {
        let d = self.diagonal();
        d.x * d.y
    }

    pub fn maximum_extent(&self) -> usize {
        let d = self.diagonal();
        if d.x > d.y {
            0
        } else {
            1
        }
    }
}

impl Index<usize> for Bounds2i {
    type Output = Point2<i32>;
    fn index(&self, i: usize) -> &Point2<i32> {
        match i {
            0 => &self.p_min,
            _ => &self.p_max,
        }
    }
}

struct Bounds2iIterator<'a> {
    point: Point2<i32>,
    bounds: &'a Bounds2i,
}

impl<'a> Iterator for Bounds2iIterator<'a> {
    type Item = &'a Point2<i32>;

    fn next(&mut self) -> Option<&'a Point2<i32>> {
        self.point.x += 1;
        if self.point.x == self.bounds.p_max.x {
            self.point.x = self.bounds.p_min.x;
            self.point.y += 1;
        }

        if self.point.y == self.bounds.p_max.y {
            None
        } else {
            Some(&self.point)
        }
    }
}

impl<'a> IntoIterator for &'a Bounds2i {
    type Item = Point2<i32>;
    type IntoIter = Bounds2iIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        Bounds2iIterator {
            point: Point2::new(self.p_min.x - 1, self.p_min.y),
            bounds: self,
        }
    }
}

impl<'a> IntoIterator for Bounds2i {
    type Item = Point2<i32>;
    type IntoIter = Bounds2iIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        Bounds2iIterator {
            point: Point2::new(self.p_min.x - 1, self.p_min.y),
            bounds: &self,
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Bounds2f {
    pub p_min: Point2<Float>,
    pub p_max: Point2<Float>,
}

impl Bounds2f {
    pub fn new(p0: Point2<Float>, p1: Point2<Float>) -> Bounds2f {
        let p_min = Point2::new(p0.x.min(p1.x), p0.y.min(p1.y));
        let p_max = Point2::new(p0.x.max(p1.x), p0.y.max(p1.y));
        Bounds2f { p_min, p_max }
    }

    pub fn corner(&self, corner: usize) -> Point2<Float> {
        let x = if corner & 1 == 0 {
            self.p_min.x
        } else {
            self.p_max.x
        };
        let y = if corner & 2 == 0 {
            self.p_min.y
        } else {
            self.p_max.y
        };

        Point2::new(x, y)
    }

    pub fn union(&self, p: &Point2<Float>) -> Bounds2f {
        let p0 = Point2::new(self.p_min.x.min(p.x), self.p_min.y.min(p.y));
        let p1 = Point2::new(self.p_max.x.max(p.x), self.p_max.y.max(p.y));
        Bounds2f {
            p_min: p0,
            p_max: p1,
        }
    }

    pub fn union_bounds(&self, b: &Bounds2f) -> Bounds2f {
        let p0 = Point2::new(self.p_min.x.min(b.p_min.x), self.p_min.y.min(b.p_min.y));
        let p1 = Point2::new(self.p_max.x.max(b.p_max.x), self.p_max.y.max(b.p_max.y));
        Bounds2f {
            p_min: p0,
            p_max: p1,
        }
    }

    pub fn intersect(&self, b: &Bounds2f) -> Bounds2f {
        let p0 = Point2::new(self.p_min.x.max(b.p_min.x), self.p_min.y.max(b.p_min.y));
        let p1 = Point2::new(self.p_max.x.min(b.p_max.x), self.p_max.y.min(b.p_max.y));
        Bounds2f {
            p_min: p0,
            p_max: p1,
        }
    }

    pub fn overlaps(&self, b: &Bounds2f) -> bool {
        (self.p_max >= b.p_min) && (self.p_min <= b.p_max)
    }

    pub fn inside(&self, p: &Point2<Float>) -> bool {
        (p >= &self.p_min) && (p <= &self.p_max)
    }

    pub fn inside_exclusive(&self, p: &Point2<Float>) -> bool {
        (p >= &self.p_min) && (p < &self.p_max)
    }

    pub fn expand(&self, delta: Float) -> Bounds2f {
        Bounds2f {
            p_min: &self.p_min - Vector2::new(delta, delta),
            p_max: &self.p_max + Vector2::new(delta, delta),
        }
    }

    pub fn diagonal(&self) -> Vector2<Float> {
        &self.p_max - &self.p_min
    }

    pub fn surface_area(&self) -> Float {
        let d = self.diagonal();
        d.x * d.y
    }

    pub fn maximum_extent(&self) -> usize {
        let d = self.diagonal();
        if d.x > d.y {
            0
        } else {
            1
        }
    }
}

impl Index<usize> for Bounds2f {
    type Output = Point2<Float>;
    fn index(&self, i: usize) -> &Point2<Float> {
        match i {
            0 => &self.p_min,
            _ => &self.p_max,
        }
    }
}

// Axis-aligned bounding box
#[derive(Debug, Copy, Clone)]
pub struct Bounds3 {
    pub p_min: Point3<Float>,
    pub p_max: Point3<Float>,
}

impl Bounds3 {
    pub fn new(p0: Point3<Float>, p1: Point3<Float>) -> Bounds3 {
        let p_min = Point3::new(p0.x.min(p1.x), p0.y.min(p1.y), p0.z.min(p1.z));
        let p_max = Point3::new(p0.x.max(p1.x), p0.y.max(p1.y), p0.z.max(p1.z));
        Bounds3 { p_min, p_max }
    }

    pub fn corner(&self, corner: usize) -> Point3<Float> {
        let x = if corner & 1 == 0 {
            self.p_min.x
        } else {
            self.p_max.x
        };
        let y = if corner & 2 == 0 {
            self.p_min.y
        } else {
            self.p_max.y
        };
        let z = if corner & 4 == 0 {
            self.p_min.z
        } else {
            self.p_max.z
        };

        Point3::new(x, y, z)
    }

    pub fn union(&self, p: &Point3<Float>) -> Bounds3 {
        let p0 = Point3::new(
            self.p_min.x.min(p.x),
            self.p_min.y.min(p.y),
            self.p_min.z.min(p.z),
        );
        let p1 = Point3::new(
            self.p_max.x.max(p.x),
            self.p_max.y.max(p.y),
            self.p_max.z.max(p.z),
        );
        Bounds3 {
            p_min: p0,
            p_max: p1,
        }
    }

    pub fn union_bounds(&self, b: &Bounds3) -> Bounds3 {
        let p0 = Point3::new(
            self.p_min.x.min(b.p_min.x),
            self.p_min.y.min(b.p_min.y),
            self.p_min.z.min(b.p_min.z),
        );
        let p1 = Point3::new(
            self.p_max.x.max(b.p_max.x),
            self.p_max.y.max(b.p_max.y),
            self.p_max.z.max(b.p_max.z),
        );
        Bounds3 {
            p_min: p0,
            p_max: p1,
        }
    }

    pub fn intersect(&self, b: &Bounds3) -> Bounds3 {
        let p0 = Point3::new(
            self.p_min.x.max(b.p_min.x),
            self.p_min.y.max(b.p_min.y),
            self.p_min.z.max(b.p_min.z),
        );
        let p1 = Point3::new(
            self.p_max.x.min(b.p_max.x),
            self.p_max.y.min(b.p_max.y),
            self.p_max.z.min(b.p_max.z),
        );
        Bounds3 {
            p_min: p0,
            p_max: p1,
        }
    }

    pub fn overlaps(&self, b: &Bounds3) -> bool {
        (self.p_max >= b.p_min) && (self.p_min <= b.p_max)
    }

    pub fn inside(&self, p: &Point3<Float>) -> bool {
        (p >= &self.p_min) && (p <= &self.p_max)
    }

    pub fn inside_exclusive(&self, p: &Point3<Float>) -> bool {
        (p >= &self.p_min) && (p < &self.p_max)
    }

    pub fn expand(&self, delta: Float) -> Bounds3 {
        Bounds3 {
            p_min: &self.p_min - Vector3::new(delta, delta, delta),
            p_max: &self.p_max + Vector3::new(delta, delta, delta),
        }
    }

    pub fn diagonal(&self) -> Vector3<Float> {
        &self.p_max - &self.p_min
    }

    pub fn area(&self) -> Float {
        let d = &self.p_max - &self.p_min;
        d.x * d.y
    }

    pub fn surface_area(&self) -> Float {
        let d = self.diagonal();
        2.0 * (d.x * d.y + d.x * d.z + d.y * d.z)
    }

    pub fn volume(&self) -> Float {
        let d = self.diagonal();
        d.x * d.y * d.z
    }

    pub fn maximum_extent(&self) -> usize {
        let d = self.diagonal();
        if d.x > d.y && d.x > d.z {
            0
        } else if d.y > d.z {
            1
        } else {
            2
        }
    }
}

impl Index<usize> for Bounds3 {
    type Output = Point3<Float>;
    fn index(&self, i: usize) -> &Point3<Float> {
        match i {
            0 => &self.p_min,
            _ => &self.p_max,
        }
    }
}
