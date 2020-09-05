use crate::bxdf_type::BxDFType;
use crate::camera::Camera;
use crate::geometry::ray::Ray;
use crate::interaction::SurfaceInteraction;
use crate::pbrt::Spectrum;
use crate::reflection::{ReflectionRecord, BSDF};
use crate::sampler::Sampler;
use crate::samplers::uniform::UniformSampler;
use crate::scene::Scene;

// Recursively evaluates radiance along reflected and refracted day directions, stopping
// at a predetermined depth.
pub struct WhittedIntegrator {
    camera: Camera,
    sampler: UniformSampler,
    max_depth: u32,
}

impl WhittedIntegrator {
    pub fn new(camera: Camera, sampler: UniformSampler, max_depth: u32) -> WhittedIntegrator {
        WhittedIntegrator {
            camera,
            sampler,
            max_depth,
        }
    }

    pub fn get_sampler(&self) -> &UniformSampler {
        &self.sampler
    }

    pub fn get_camera(&self) -> &Camera {
        &self.camera
    }

    // Returns the radiance arriving at the origin of the given ray
    //TODO make generic sampler
    pub fn li(
        &self,
        ray: &Ray,
        scene: &Scene,
        sampler: &mut UniformSampler,
        depth: u32,
    ) -> Spectrum {
        // Find closes ray intersection or return background radiance
        match scene.intersect(ray) {
            Some(mut isect) => {
                // compute emitted and reflected light at ray intersection point
                let normal = isect.interaction.normal;
                let w_out = isect.interaction.w0;

                isect.compute_scattering_functions(ray);

                // compute emitted light ray if ray hit an rea light source
                let mut lighting = isect.le(&w_out);

                // TODO BSDF is not used!
                if let Some(scattering_function) = isect
                    // .shape
                    // .get_material()
                    .material
                    .compute_scattering_functions(&isect, false)
                {
                    // add contribution of each light source
                    for light in scene.get_lights() {
                        let scatter_record = light.sample_li(&isect, &sampler.get_2d());
                        if scatter_record.illuminance.is_black() || scatter_record.pdf == 0.0 {
                            continue;
                        }
                        let f =
                            scattering_function.f(&w_out, &scatter_record.wi, BxDFType::BsdfAll);
                        if !f.is_black() && scatter_record.visibility_tester.unoccluded(scene) {
                            lighting += f
                                * scatter_record.illuminance
                                * scatter_record.wi.dot(&normal).abs()
                                / scatter_record.pdf;
                        }
                    }

                    if depth + 1 < self.max_depth {
                        // trace rays for specular reflection and refraction
                        lighting += self.specular_reflect(
                            &ray,
                            &isect,
                            &scattering_function,
                            scene,
                            sampler,
                            depth,
                        );
                        lighting += self.specular_transmit(
                            &ray,
                            &isect,
                            &scattering_function,
                            scene,
                            sampler,
                            depth,
                        );
                    }
                }

                return lighting;
            }
            _ => {
                let mut spectrum = Spectrum::new(0.0);
                for light in scene.get_lights() {
                    spectrum += light.le(ray);
                }
                return spectrum;
            }
        }
    }

    fn specular_reflect(
        &self,
        ray: &Ray,
        isect: &SurfaceInteraction,
        bsdf: &BSDF,
        scene: &Scene,
        sampler: &mut UniformSampler,
        depth: u32,
    ) -> Spectrum {
        // compute specular reflection direction wi and BSDF value
        let wo = &isect.interaction.w0;
        let bsdf_flags = BxDFType::BsdfReflection | BxDFType::BsdfSpecular;
        let mut sample_type = BxDFType::Empty; //TODO
        let reflection_record = bsdf.sample_f(wo, &sampler.get_2d(), bsdf_flags, sample_type);

        // return contribution of specular reflection
        let ns = &isect.shading.normal;
        match reflection_record {
            Some(reflection_record) => {
                if reflection_record.pdf > 0.0
                    && !reflection_record.illuminance.is_black()
                    && reflection_record.wi.dot(ns).abs() != 0.0
                {
                    // compute ray differential rd for specular reflection
                    let mut rd = isect.interaction.spawn_ray(&reflection_record.wi);
                    rd = match &ray.differential {
                        Some(diff) => {
                            let rx_origin = &isect.interaction.point + &isect.dpdx;
                            let ry_origin = &isect.interaction.point + &isect.dpdy;
                            // compute differential reflected directions
                            let dndx = &isect.shading.dndu * isect.dudx
                                + (&isect.shading.dndv * isect.dvdx);
                            let dndy =
                                &isect.shading.dndu * isect.dudy + &isect.shading.dndv * isect.dvdy;
                            let dwodx = -&diff.rx_direction - wo;
                            let dwody = -&diff.ry_direction - wo;
                            let dDNdx = dwodx.dot(&ns) + wo.dot(&dndx);
                            let dDNdy = dwody.dot(&ns) + wo.dot(&dndy);
                            let rx_direction = &reflection_record.wi - dwodx
                                + 2.0 * (wo.dot(&ns) * dndx + dDNdx * ns);
                            let ry_direction = &reflection_record.wi - dwody
                                + 2.0 * (wo.dot(&ns) * dndy + dDNdy * ns);

                            rd.with_differential(rx_origin, ry_origin, rx_direction, ry_direction)
                        }
                        _ => rd,
                    };
                    // TODO can we just return the reflection record and make the recursive call in
                    // the caller function
                    reflection_record.illuminance
                        * self.li(&rd, scene, sampler, depth + 1)
                        * reflection_record.wi.dot(&ns)
                        / reflection_record.pdf
                } else {
                    Spectrum::new(0.0)
                }
            }
            None => Spectrum::new(0.0),
        }
    }

    //Requests transmision specular component of bsdf rather than reflection component
    fn specular_transmit(
        &self,
        ray: &Ray,
        isect: &SurfaceInteraction,
        bsdf: &BSDF,
        scene: &Scene,
        sampler: &mut UniformSampler,
        depth: u32,
    ) -> Spectrum {
        // compute specular reflection direction wi and BSDF value
        let wo = &isect.interaction.w0;
        let bsdf_flags = BxDFType::BsdfTransmission | BxDFType::BsdfSpecular;
        let mut sample_type = BxDFType::Empty; //TODO
        let reflection_record = bsdf.sample_f(wo, &sampler.get_2d(), bsdf_flags, sample_type);

        // return contribution of specular reflection
        let ns = &isect.shading.normal;
        match reflection_record {
            Some(reflection_record) => {
                if reflection_record.pdf > 0.0
                    && !reflection_record.illuminance.is_black()
                    && reflection_record.wi.dot(ns).abs() != 0.0
                {
                    // compute ray differential rd for specular reflection
                    let mut rd = isect.interaction.spawn_ray(&reflection_record.wi);
                    rd = match &ray.differential {
                        Some(diff) => {
                            let rx_origin = &isect.interaction.point + &isect.dpdx;
                            let ry_origin = &isect.interaction.point + &isect.dpdy;
                            // compute differential reflected directions
                            let dndx =
                                &isect.shading.dndu * isect.dudx + &isect.shading.dndv * isect.dvdx;
                            let dndy =
                                &isect.shading.dndu * isect.dudy + &isect.shading.dndv * isect.dvdy;
                            let dwodx = -&diff.rx_direction - wo;
                            let dwody = -&diff.ry_direction - wo;
                            let dDNdx = dwodx.dot(&ns) + wo.dot(&dndx);
                            let dDNdy = dwody.dot(&ns) + wo.dot(&dndy);
                            let rx_direction = &reflection_record.wi - dwodx
                                + 2.0 * (wo.dot(&ns) * dndx + dDNdx * ns);
                            let ry_direction = &reflection_record.wi - dwody
                                + 2.0 * (wo.dot(&ns) * dndy + dDNdy * ns);

                            rd.with_differential(rx_origin, ry_origin, rx_direction, ry_direction)
                        }
                        _ => rd,
                    };
                    // TODO can we just return the reflection record and make the recursive call in
                    // the caller function
                    reflection_record.illuminance
                        * self.li(&rd, scene, sampler, depth + 1)
                        * reflection_record.wi.dot(&ns)
                        / reflection_record.pdf
                } else {
                    Spectrum::new(0.0)
                }
            }
            None => Spectrum::new(0.0),
        }
    }
}
