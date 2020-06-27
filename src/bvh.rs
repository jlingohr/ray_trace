use super::aabb::AABB;
use super::hittable::{HitRecord, Hittable};
use super::ray::Ray;

use rand::Rng;
use std::cmp::Ordering;

pub struct BVH {
    tree: BVHNode,
    bbox: AABB,
}

enum BVHNode {
    Branch { left: Box<BVH>, right: Box<BVH> },
    Leaf(Box<dyn Hittable>),
}

impl BVH {
    pub fn new(mut objects: Vec<Box<dyn Hittable>>, time0: f64, time1: f64) -> BVH {
        fn box_compare(
            time0: f64,
            time1: f64,
            axis: usize,
        ) -> impl FnMut(&Box<dyn Hittable>, &Box<dyn Hittable>) -> Ordering {
            move |a, b| {
                let box_a = a.bounding_box(time0, time1);
                let box_b = b.bounding_box(time0, time1);
                if let (Some(a), Some(b)) = (box_a, box_b) {
                    let ac = a.min()[axis] + a.max()[axis];
                    let bc = b.min()[axis] + b.max()[axis];
                    ac.partial_cmp(&bc).unwrap()
                } else {
                    panic!("No bounding box in bvh node")
                }
            }
        }

        fn axis_range(
            hittable: &Vec<Box<dyn Hittable>>,
            time0: f64,
            time1: f64,
            axis: usize,
        ) -> f64 {
            let (min, max) = hittable
                .iter()
                .fold((f64::MAX, f64::MIN), |(bmin, bmax), hit| {
                    if let Some(aabb) = hit.bounding_box(time0, time1) {
                        (bmin.min(aabb.min()[axis]), bmax.max(aabb.max()[axis]))
                    } else {
                        (bmin, bmax)
                    }
                });
            max - min
        }

        let mut axis_ranges: Vec<(usize, f64)> = (0..3)
            .map(|a| (a, axis_range(&objects, time0, time1, a)))
            .collect();

        axis_ranges.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        let axis = axis_ranges[0].0;
        objects.sort_unstable_by(box_compare(time0, time1, axis));
        let len = objects.len();

        match len {
            0 => panic!("no elements in scene"),
            1 => {
                let leaf = objects.pop().unwrap();
                if let Some(bbox) = leaf.bounding_box(time0, time1) {
                    BVH {
                        tree: BVHNode::Leaf(leaf),
                        bbox,
                    }
                } else {
                    panic!("No bounding box in bvh node")
                }
            }
            _ => {
                let right = BVH::new(objects.drain(len / 2..).collect(), time0, time1);
                let left = BVH::new(objects, time0, time1);
                let bbox = AABB::surrounding_box(&left.bbox, &right.bbox);
                BVH {
                    tree: BVHNode::Branch {
                        left: Box::new(left),
                        right: Box::new(right),
                    },
                    bbox,
                }
            }
        }
    }
}

impl Hittable for BVH {
    fn hit(&self, ray: &Ray, t_min: f64, mut t_max: f64) -> Option<HitRecord> {
        if self.bbox.hit(ray, t_min, t_max) {
            match &self.tree {
                BVHNode::Leaf(leaf) => leaf.hit(ray, t_min, t_max),
                BVHNode::Branch { left, right } => {
                    let left = left.hit(ray, t_min, t_max);
                    if let Some(l) = &left {
                        t_max = l.t
                    };
                    let right = right.hit(ray, t_min, t_max);
                    if right.is_some() {
                        right
                    } else {
                        left
                    }
                }
            }
        } else {
            None
        }
    }

    fn bounding_box(&self, _t0: f64, _t1: f64) -> Option<AABB> {
        Some(self.bbox)
    }
}
