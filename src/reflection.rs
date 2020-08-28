use crate::pbrt::{clamp, Float, Spectrum, INV_PI};
use nalgebra::{Point2, Vector3};
use std::cmp::max;

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
trait BxDF {
    fn matches_flags(&self, t: BxDFType) -> bool {
        self.get_type() & t == self.get_type()
    }

    // Returns the valie of the distribution function for the given pair of directions.
    // Implicity assumes that light in different wavelengths is decoupled
    fn f(&self, wo: &Vector3<Float>, wi: &Vector3<Float>) -> Spectrum;

    // Handle scattering described by delta distributions as well as for randomly sampling directions
    // that scatter light along multiple directions.
    fn sample_f(
        &self,
        wo: &Vector3<Float>,
        wi: &Vector3<Float>,
        sample: &Point2<Float>,
        pdf: Float,
    ) -> Spectrum;

    // computes hemispherical-directional reflectance function, either in closed form or using Monte Carlo integration
    fn rho_hdr(&self, wo: &Vector3<Float>, n_samples: u32, samples: &Point2<Float>) -> Spectrum;

    // compute the hemispherical-hemispherical reflectance of a surface
    fn rho_hhr(
        &self,
        n_samples: u32,
        samples1: &Point2<Float>,
        samples2: &Point2<Float>,
    ) -> Spectrum;

    fn get_type(&self) -> BxDFType;
}

// odels a perfect diffuse surface that scatters incident illumination equally in all directions
pub struct LambertianReflection {
    reflectance: Spectrum,
}

impl LambertianReflection {
    pub fn new(reflectance: Spectrum) -> LambertianReflection {
        LambertianReflection { reflectance }
    }
}

impl BxDF for LambertianReflection {
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

    fn rho_hdr(&self, wo: &Vector3<f64>, n_samples: u32, samples: &Point2<f64>) -> Spectrum {
        self.reflectance
    }

    fn rho_hhr(&self, n_samples: u32, samples1: &Point2<f64>, samples2: &Point2<f64>) -> Spectrum {
        self.reflectance
    }

    fn get_type(&self) -> BxDFType {
        BxDFType::BsdfReflection | BxDFType::BsdfDiffuse
    }
}
