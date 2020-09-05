use crate::interaction::SurfaceInteraction;
use crate::pbrt::{clamp, Float, Spectrum};
use crate::reflection::{BSDF, BxDF, LambertianReflection};
use crate::texture::Texture;

#[derive(Debug, Copy, Clone)]
pub enum Material {
    Matte(MatteMaterial),
}

impl Material {
    // Given object containing geometric properties at an intersection point on the surface,
    // determines the reflective properties at the point and initializes bsdf member variables
    pub fn compute_scattering_functions(
        &self,
        isect: &SurfaceInteraction,
        allow_multiple_lobes: bool,
    ) -> Option<BSDF> {
        match self {
            Material::Matte(material) => material.compute_scattering_functions(isect),
        }
    }

    //
    pub fn bump(&self, d: &Texture<Float>, isect: &SurfaceInteraction) {
        unimplemented!()
    }
}

// MatteMaterial describes a purely diffuse surface
// It is parameterized by a spectral diffuse reflection value and a scalar roughness value.
// Takes optional scalat texture bump_map that defines an offset function over the surface.
// If not None, this texture is used to compute a shadig normal at each point based on
// the function it defines.
#[derive(Debug, Copy, Clone)]
pub struct MatteMaterial {
    kd: Texture<Spectrum>,
    sigma: Texture<Float>,
    bump_map: Option<Texture<Float>>,
}

impl MatteMaterial {
    pub fn new(
        kd: Texture<Spectrum>,
        sigma: Texture<Float>,
        bump_map: Option<Texture<Float>>,
    ) -> MatteMaterial {
        MatteMaterial {
            kd,
            sigma,
            bump_map,
        }
    }

    // Determine bump map's effect on the shading geometry, evaluates textures, and
    // allocates and returning the appropriate bsdf
    fn compute_scattering_functions(&self, isect: &SurfaceInteraction) -> Option<BSDF> {
        // if bump_map present, calculate the shading normal at the point
        //TODO add this
        // self.bump_map.map(|bump_map| bump(bump_map, isect));

        // evaluate textures for MatteMaterial and allocate brdf
        let mut bsdf = BSDF::new(isect, 1.0);
        let r = self.kd.evaluate(isect).clamp(0.0, 1.0);
        let sig = clamp(*self.sigma.evaluate(isect), 0.0, 90.0);
        if r.is_black() {
            if sig == 0.0 {
                bsdf.add(BxDF::LambertianRef(LambertianReflection::new(r)));
            } else {
                bsdf.add(BxDF::LambertianRef(LambertianReflection::new(r))); //TODO just for now
                // bsdf.add(oren_nayar(r, sig));
            }
        }
        Some(bsdf)
    }

    fn bump(&self, d: &Texture<Float>, isect: &SurfaceInteraction) {
        unimplemented!()
    }
}
