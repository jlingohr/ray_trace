use std::sync::Arc;

use crate::geometry::bounds::Bounds3;
use crate::geometry::ray::Ray;
use crate::interaction::SurfaceInteraction;
use crate::light::Light;
use crate::primitive::Aggregate;

pub struct Scene {
    lights: Vec<Light>,
    aggregate: Arc<Aggregate>,
    world_bound: Bounds3,
}

impl Scene {
    pub fn new(aggregate: Arc<Aggregate>, lights: Vec<Light>) -> Scene {
        //TODO need to preprocess lights?
        let world_bound = aggregate.world_bound();
        Scene {
            lights,
            aggregate,
            world_bound,
        }
    }

    pub fn get_lights(&self) -> &Vec<Light> {
        &self.lights
    }

    pub fn world_bound(&self) -> Bounds3 {
        unimplemented!()
    }

    // Traces the given ray into the scene and returns nearest Some<SurfaceInteraction> if the ray
    // intersected any of the primitives, otherwise returns None
    pub fn intersect(&self, ray: &Ray) -> Option<SurfaceInteraction> {
        self.aggregate.intersect(ray)
    }

    // Checks for the existence of intersections along a ray
    pub fn intersect_p(&self, ray: &Ray) -> bool {
        self.aggregate.intersect_p(ray)
    }
}
