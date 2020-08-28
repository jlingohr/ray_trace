use crate::interaction::SurfaceInteraction;
use crate::pbrt::{clamp, Float, Spectrum, INV_PI};
use crate::reflection::BxDFType::{BsdfReflection, BsdfTransmission};
use nalgebra::{Point2, Vector3};
use std::cmp::max;

const MAX_BXDFS: usize = 8;

// Uitility functions to compute angles for the shading coordinate system

#[inline]
fn cos_theta(w: &Vector3<Float>) -> Float {
    w.z
}

#[inline]
fn cos_2_theta(w: &Vector3<Float>) -> Float {
    w.z * w.z
}

#[inline]
fn abs_cos_theta(w: &Vector3<Float>) -> Float {
    w.z.abs()
}

#[inline]
fn sin_2_theta(w: &Vector3<Float>) -> Float {
    max(0.0, 1.0 - cos_2_theta(w))
}

#[inline]
fn sin_theta(w: &Vector3<Float>) -> Float {
    sin_2_theta(w).sqrt()
}

#[inline]
fn tan_theta(w: &Vector3<Float>) -> Float {
    sin_theta(w) / cos_theta(w)
}

#[inline]
fn tan_2_theta(w: &Vector3<Float>) -> Float {
    sin_2_theta(w) / cos_theta(w)
}

#[inline]
fn cos_phi(w: &Vector3<Float>) -> Float {
    let sin_theta = sin_theta(w);
    if sin_theta == 0.0 {
        1.0
    } else {
        clamp(w.x / sin_theta, -1.0, 1.0)
    }
}

#[inline]
fn sin_phi(w: &Vector3<Float>) -> Float {
    let sin_theta = sin_theta(w);
    if sin_theta == 0.0 {
        1.0
    } else {
        clamp(w.y / sin_theta, -1.0, 1.0)
    }
}

#[inline]
fn cos_2_phi(w: &Vector3<Float>) -> Float {
    cos_phi(w) * cos_phi(w)
}

#[inline]
fn sin_2_phi(w: &Vector3<Float>) -> Float {
    sin_phi(w) * sin_phi(w)
}

#[inline]
fn cos_d_phi(wa: &Vector3<Float>, wb: &Vector3<Float>) -> Float {
    clamp(
        (wa.x * wb.x + wa.y * wb.y) / (wa.x * wa.x + wa.y * wa.y).sqrt()
            * (wb.x * wb.x + wb.y * wb.y),
        -1.0,
        1.0,
    )
}

// For each BxDF the flags should have at least one of REFLECTION or TRANSMISSION set
// and exactly one of diffuse, glossy, and specular flags.
#[repr(u8)]
#[derive(Clone, Copy)]
enum BxDFType {
    BsdfReflection = 1,
    BsdfTransmission = 2,
    BsdfDiffuse = 4,
    BsdfGlossy = 8,
    BsdfSpecular = 16,
    BsdfAll = 31,
}

//TODO not all BxDF can use f or sample_f so probably refactor these out into another trait
// and a function called pdf
pub enum BxDF {
    Empty,
    LambertianRef(LambertianReflection),
}

impl BxDF {
    fn matches_flags(&self, t: BxDFType) -> bool {
        self.get_type() & t == self.get_type()
    }

    // Returns the valie of the distribution function for the given pair of directions.
    // Implicity assumes that light in different wavelengths is decoupled
    fn f(&self, wo: &Vector3<Float>, wi: &Vector3<Float>) -> Spectrum {
        unimplemented!()
    }

    // Handle scattering described by delta distributions as well as for randomly sampling directions
    // that scatter light along multiple directions.
    fn sample_f(
        &self,
        wo: &Vector3<Float>,
        wi: &Vector3<Float>,
        sample: &Point2<Float>,
        pdf: Float,
    ) -> Spectrum {
        unimplemented!()
    }

    // computes hemispherical-directional reflectance function, either in closed form or using Monte Carlo integration
    fn rho_hdr(
        &self,
        wo: &Vector3<Float>,
        n_samples: usize,
        samples: &Vec<Point2<Float>>,
    ) -> Spectrum {
        unimplemented!()
    }

    // compute the hemispherical-hemispherical reflectance of a surface
    fn rho_hhr(
        &self,
        n_samples: usize,
        samples1: &Vec<Point2<Float>>,
        samples2: &Vec<Point2<Float>>,
    ) -> Spectrum {
        unimplemented!()
    }

    fn get_type(&self) -> BxDFType {
        unimplemented!()
    }
}

// models a perfect diffuse surface that scatters incident illumination equally in all directions
pub struct LambertianReflection {
    reflectance: Spectrum,
}

impl LambertianReflection {
    pub fn new(reflectance: Spectrum) -> LambertianReflection {
        LambertianReflection { reflectance }
    }

    fn f(&self, wo: &Vector3<f64>, wi: &Vector3<f64>) -> Spectrum {
        INV_PI * self.reflectance
    }

    fn sample_f(
        &self,
        wo: &Vector3<f64>,
        wi: &Vector3<f64>,
        sample: &Point2<f64>,
        pdf: f64,
    ) -> Spectrum {
        unimplemented!()
    }

    fn rho_hdr(&self, wo: &Vector3<Float>, n_samples: usize, samples: &Point2<Float>) -> Spectrum {
        self.reflectance
    }

    fn rho_hhr(
        &self,
        n_samples: u32,
        samples1: Vec<Point2<Float>>,
        samples2: &Point2<Float>,
    ) -> Spectrum {
        self.reflectance
    }

    fn get_type(&self) -> BxDFType {
        BxDFType::BsdfReflection | BxDFType::BsdfDiffuse
    }
}

// Represent a collection of BRDFs and BTDFs
pub struct BSDF {
    pub eta: Float,
    normal_shading: Vector3<Float>,
    normal_geometric: Vector3<Float>,
    ss: Vector3<Float>,
    ts: Vector3<Float>,
    num_bxdfs: usize,
    bxdfs: [BxDF; MAX_BXDFS],
}

impl BSDF {
    // SurfaceInteraction contains information about the differential geometry at the point
    // on a surface. Eta gives the relative index of refraction over the boundary.
    // For surfaces where eta isn't used it should be 1.
    pub fn new(isect: &SurfaceInteraction, eta: Float) -> BSDF {
        let normal_shading = isect.shading.normal.clone_owned();
        let normal_geometric = isect.interaction.normal.clone_owned();
        let ss = isect.shading.dpdu.normalize();
        let ts = normal_geometric.cross(&ss);

        BSDF {
            eta,
            normal_shading,
            normal_geometric,
            ss,
            ts,
            num_bxdfs: 0,
            bxdfs: [BxDF::Empty; MAX_BXDFS],
        }
    }

    // Return numner of BxDFs stored by the BSDF that match a particular set of flags
    pub fn num_components(&self, flags: BxDFType) -> i32 {
        let mut count = 0;
        for i in 0..self.num_bxdfs {
            if self.bxdfs[i].matches_flags(flags) {
                count += 1
            }
        }
        count
    }

    // Perform transforms to local coordinate system
    pub fn world_to_local(&self, v: &Vector3<Float>) -> Vector3<Float> {
        Vector3::new(
            v.dot(&self.ss),
            v.dot(&self.ts),
            v.dot(&self.normal_shading),
        )
    }

    pub fn local_to_world(&self, v: &Vector3<Float>) -> Vector3<Float> {
        Vector3::new(
            self.ss.x * v.x + self.ts.x * v.y + self.normal_shading.x * v.z,
            self.ss.y * v.x + self.ts.y * v.y + self.normal_shading.y * v.z,
            self.ss.z * v.x + self.ts.z * v.y + self.normal_shading.z * v.z,
        )
    }

    pub fn f(&self, wo_w: &Vector3<Float>, wi_w: &Vector3<Float>, flags: BxDFType) -> Spectrum {
        let wi = self.world_to_local(wi_w);
        let wo = self.world_to_local(wo_w);
        let reflect = wi_w.dot(&self.normal_geometric) & wo_w.dot(&self.normal_geometric) > 0.0;

        let mut spectrum = Spectrum::new(0.0);
        for i in 0..self.num_bxdfs {
            if self.bxdfs[i].matches_flags(flags)
                && ((reflect && (self.bxdfs[i].get_type() * BsdfReflection))
                    || (!reflect && (self.bxdfs[i].get_type() & BsdfTransmission)))
            {
                spectrum += self.bxdfs[i].f(&wo, &wi);
            }
        }
        spectrum
    }

    pub fn rho_hdr(
        &self,
        wo: &Vector3<Float>,
        n_samples: usize,
        samples: &Vec<Point2<Float>>,
        flags: BxDFType,
    ) -> Spectrum {
        let world_to_local = self.world_to_local(wo);
        let mut spectrum = Spectrum::new(0.0);
        for i in 0..self.num_bxdfs {
            if self.bxdfs[i].matches_flags(flags) {
                spectrum += self.bxdfs[i].rho_hdr(wo, n_samples, &samples);
            }
        }
        spectrum
    }

    pub fn rho_hhr(
        &self,
        n_samples: usize,
        samples1: &Vec<Point2<Float>>,
        samples2: &Vec<Point2<Float>>,
        flags: BxDFType,
    ) -> Spectrum {
        let mut spectrum = Spectrum::new(0.0);
        for i in 0..self.num_bxdfs {
            if self.bxdfs[i].matches_flags(flags) {
                spectrum += self.bxdfs[i].rho_hhr(n_samples, samples1, samples2);
            }
        }
        spectrum
    }

    pub fn add(&mut self, b: BxDF) {
        if self.num_bxdfs == MAX_BXDFS {
            //TODO should handle this correctly
            panic!("Maximum number of BxDFs reached")
        }
        self.bxdfs[self.num_bxdfs] = b;
        self.num_bxdfs += 1;
    }
}
