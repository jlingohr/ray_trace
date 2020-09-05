use std::cmp::max;

use nalgebra::{distance_squared, Point2, Point3, Vector3};

use crate::geometry::ray::Ray;
use crate::interaction::{CommonInteraction, SurfaceInteraction};
use crate::pbrt::{Float, Spectrum};
use crate::samplers::uniform::UniformSampler;
use crate::scene::Scene;
use crate::transforms::Transform;

// Struct to store light scattering record
pub struct ScatterRecord {
    pub wi: Vector3<Float>,
    pub pdf: Float,
    pub visibility_tester: VisibilityTester,
    pub illuminance: Spectrum,
}

// Represent flags for flags mask field characterizing varius kinds of light sources
#[repr(u32)]
#[derive(Debug, Copy, Clone)]
pub enum LightFlags {
    DeltaPosition = 1,
    DeltaDirection = 2,
    Area = 4,
    Infinite = 8,
}

impl LightFlags {
    pub fn is_delta_flag(self) -> bool {
        match self {
            LightFlags::DeltaPosition => true,
            LightFlags::DeltaDirection => true,
            _ => false,
        }
    }
}

// Closure that allows lights to return a radiance value under the assumption
// that the reference point and the light source are mutually visible. This allows the
// integrator to decide if illumination from the incident direction is relevant before
// incurring the cost of tracing the shadow ray.
//TODO how to use actual Rust closure?
#[derive(Copy, Clone)]
pub struct VisibilityTester {
    p0: CommonInteraction,
    p1: CommonInteraction,
}

impl VisibilityTester {
    pub fn new(p0: CommonInteraction, p1: CommonInteraction) -> VisibilityTester {
        VisibilityTester { p0, p1 }
    }

    // Traces a shadow ray between two points and returns bool result, ignoring the effects
    // of any scattering medium that the ray passes through
    pub fn unoccluded(&self, scene: &Scene) -> bool {
        scene.intersect_p(&self.p0.spawn_ray_from_interaction(&self.p1))
    }

    // Compute beam transmittance when we don't want to ignore effects of scattering medium
    // TODO should be Sampler enum or trait
    pub fn transmittance(&self, scene: &Scene, sampler: &UniformSampler) -> Spectrum {
        let mut ray = self.p0.spawn_ray_from_interaction(&self.p1);
        let mut transmittance = Spectrum::new(1.0);

        loop {
            if let Some(hit_surface) = scene.intersect(&ray) {
                // if hit_surface.material {
                //     Spectrum::new(0.0)
                // }
                //TODO
                // if ray.medium {
                //     transmittance += medium.trace(&ray, sampler);
                // }
                ray = hit_surface.interaction.spawn_ray_from_interaction(&self.p1);
            } else {
                break;
            }
        }
        transmittance
    }
}

pub enum Light {
    PointLight(PointLight),
}

impl Light {
    pub fn sample_li(
        &self,
        interaction: &SurfaceInteraction, //TODO need trait later becaue of MediumInteraction
        u: &Point2<Float>,
    ) -> ScatterRecord {
        unimplemented!()
    }

    pub fn le(&self, ray: &Ray) -> Spectrum {
        unimplemented!()
    }

    // Given an interaction providing the world space position of a reference point in a scene and
    // a time associated with it, return the radiance at that time due to that light
    pub fn pdf_li(&self, interaction: &SurfaceInteraction, wi: &Vector3<Float>) {
        unimplemented!()
    }

    pub fn sample_le(
        &self,
        ul: &Point2<Float>,
        u2: &Point2<Float>,
        time: Float,
        ray: &Ray,
        nlight: &Vector3<Float>,
        pdf_pos: Float,
        pdf_df: Float,
    ) -> Spectrum {
        unimplemented!()
    }
    fn pdf_le(&self, ray: &Ray, nlight: &Vector3<Float>, pdf_pos: Float, pdf_dir: Float) {
        unimplemented!()
    }
}

pub struct BaseLight {
    pub flags: LightFlags,
    pub n_samples: i32,
    // pub medium_interface: MediumInterface
    pub light_to_world: Transform,
    pub world_to_light: Transform,
}

impl BaseLight {
    pub fn new(flags: LightFlags, light_to_world: Transform, n_samples: i32) -> BaseLight {
        let n_samples = max(1, n_samples);
        let world_to_light = light_to_world.inverse();
        BaseLight {
            flags,
            n_samples,
            light_to_world,
            world_to_light,
        }
    }
}

pub struct AreaLight {
    base_light: BaseLight,
}

impl AreaLight {
    pub fn new(light_to_world: Transform, n_samples: i32) -> AreaLight {
        let base_light = BaseLight::new(LightFlags::Area, light_to_world, n_samples);
        AreaLight { base_light }
    }
}

// PointLight represents an isotropic point light source that emits the same amount of light
// in all directions.
pub struct PointLight {
    base_light: BaseLight,
    position: Point3<Float>,
    intensity: Spectrum,
}

impl PointLight {
    pub fn new(light_to_world: Transform, intensity: Spectrum) -> PointLight {
        let base_light = BaseLight::new(LightFlags::DeltaPosition, light_to_world, 1);
        let position = light_to_world.transform_point(&Point3::new(0.0, 0.0, 0.0));

        PointLight {
            base_light,
            position,
            intensity,
        }
    }

    // Report illumination arriving at a point for all types of light sources
    pub fn sample_li(&self, interaction: CommonInteraction, u: &Point2<Float>) -> ScatterRecord {
        let wi = (&self.position - &interaction.point).normalize();
        let pdf = 1.0;
        let common_interaction = CommonInteraction::new(
            self.position.clone(),
            interaction.time,
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 0.0),
        );
        let visibility_tester = VisibilityTester::new(interaction, common_interaction);

        let illuminance = self.intensity / distance_squared(&self.position, &interaction.point);

        ScatterRecord {
            wi,
            pdf,
            visibility_tester,
            illuminance,
        }
    }
}
