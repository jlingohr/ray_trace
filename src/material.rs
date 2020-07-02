use super::color::Color;
use super::hittable::HitRecord;
use super::pdf;
use super::pdf::PDF;
use super::ray::Ray;
use super::texture::Texture;
use super::utils;
use super::vector::Vec3;
use rand::prelude::ThreadRng;

pub enum ScatterRecord<'a> {
    Specular { scattered_ray: Ray, albedo: Color },
    Scatter { albedo: Color, pdf: PDF<'a> },
}

pub trait Material: Sync {
    fn scatter(&self, ray: Ray, rec: &HitRecord, rng: &mut ThreadRng) -> Option<ScatterRecord>;
    fn emitted(&self, _ray: Ray, _rec: &HitRecord) -> Color {
        Color::new(0.0, 0.0, 0.0)
    }
    fn scattering_pdf(&self, ray: &Ray, rec: &HitRecord, scattered: &Ray) -> f64 {
        1.0
    }
}

#[derive(Copy, Clone)]
pub struct Lambertian<T: Texture> {
    pub albedo: T,
}

impl<T: Texture> Lambertian<T> {
    pub fn new(albedo: T) -> Lambertian<T> {
        Lambertian { albedo }
    }
}

impl<T: Texture> Material for Lambertian<T> {
    fn scatter(&self, ray: Ray, rec: &HitRecord, rng: &mut ThreadRng) -> Option<ScatterRecord> {
        Some(ScatterRecord::Scatter {
            albedo: self.albedo.value(rec.u, rec.v, &rec.point),
            pdf: PDF::cosine_pdf(&rec.normal),
        })
    }

    fn scattering_pdf(&self, ray: &Ray, rec: &HitRecord, scattered: &Ray) -> f64 {
        let cosine = rec.normal.dot(&scattered.direction.unit()).max(0.0);
        cosine / std::f64::consts::PI
    }
}

#[derive(Copy, Clone)]
pub struct Metal {
    pub albedo: Color,
    fuzz: f64,
}

impl Metal {
    pub fn new(albedo: Color, fuzz: f64) -> Metal {
        Metal { albedo, fuzz }
    }
}

impl Material for Metal {
    fn scatter(&self, ray: Ray, rec: &HitRecord, rng: &mut ThreadRng) -> Option<ScatterRecord> {
        let reflected = ray.direction.unit().reflect(rec.normal);
        let scattered = Ray::new(
            rec.point,
            reflected + self.fuzz * Vec3::random_in_unit_sphere(rng),
            ray.time,
        );
        if scattered.direction.dot(&rec.normal) > 0.0 {
            Some(ScatterRecord::Specular {
                scattered_ray: scattered,
                albedo: self.albedo,
            })
        } else {
            None
        }
    }
}

#[derive(Copy, Clone)]
pub struct Dielectric {
    pub ref_idx: f64,
}

impl Dielectric {
    pub fn new(ref_idx: f64) -> Dielectric {
        Dielectric { ref_idx }
    }

    fn schlick(cosine: f64, ref_idx: f64) -> f64 {
        let r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
        let r0 = r0 * r0;
        r0 + (1.0 - r0) * (1.0 - cosine).powi(5)
    }
}

impl Material for Dielectric {
    fn scatter(&self, ray: Ray, rec: &HitRecord, rng: &mut ThreadRng) -> Option<ScatterRecord> {
        let attenuation = Color::new(1.0, 1.0, 1.0);
        let etai_over_etat = if rec.front_face {
            1.0 / self.ref_idx
        } else {
            self.ref_idx
        };
        let unit_direction = ray.direction.unit();
        let cos_theta = f64::min(-unit_direction.dot(&rec.normal), 1.0);
        let sin_theta = (1.0 - (cos_theta * cos_theta)).sqrt();
        let scattered = if etai_over_etat * sin_theta > 1.0 {
            let reflected = unit_direction.reflect(rec.normal);
            ScatterRecord::Specular {
                scattered_ray: Ray::new(rec.point, reflected, ray.time),
                albedo: attenuation,
            }
        } else if utils::random_double(rng) < Dielectric::schlick(cos_theta, etai_over_etat) {
            let reflected = unit_direction.reflect(rec.normal);
            ScatterRecord::Specular {
                scattered_ray: Ray::new(rec.point, reflected, ray.time),
                albedo: attenuation,
            }
        } else {
            let refracted = unit_direction.refract(rec.normal, etai_over_etat);
            ScatterRecord::Specular {
                scattered_ray: Ray::new(rec.point, refracted, ray.time),
                albedo: attenuation,
            }
        };
        Some(scattered)
    }
}

#[derive(Copy, Clone)]
pub struct DiffuseLight<T: Texture> {
    pub emit: T,
}

impl<T: Texture> DiffuseLight<T> {
    pub fn new(texture: T) -> DiffuseLight<T> {
        DiffuseLight { emit: texture }
    }
}

impl<T: Texture> Material for DiffuseLight<T> {
    fn scatter(&self, _ray: Ray, _rec: &HitRecord, _rng: &mut ThreadRng) -> Option<ScatterRecord> {
        None
    }

    fn emitted(&self, ray: Ray, rec: &HitRecord) -> Color {
        if rec.front_face {
            self.emit.value(rec.u, rec.v, &rec.point)
        } else {
            Color::new(0.0, 0.0, 0.0)
        }
    }
}

#[derive(Copy, Clone)]
pub struct Isotropic<T: Texture> {
    albedo: T,
}

impl<T: Texture> Isotropic<T> {
    pub fn new(a: T) -> Isotropic<T> {
        Isotropic { albedo: a }
    }
}

impl<T: Texture> Material for Isotropic<T> {
    fn scatter(&self, ray: Ray, rec: &HitRecord, rng: &mut ThreadRng) -> Option<ScatterRecord> {
        let scattered = Ray::new(rec.point, Vec3::random_in_unit_sphere(rng), ray.time);
        let attenuation = self.albedo.value(rec.u, rec.v, &rec.point);
        Some(ScatterRecord::Specular {
            scattered_ray: scattered,
            albedo: attenuation,
        })
    }
}
