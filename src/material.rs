use super::color::Color;
use super::hittable::HitRecord;
use super::ray::Ray;
use super::texture::Texture;
use super::utils;
use super::vector::Vec3;
use rand::prelude::ThreadRng;

pub struct Scatter {
    pub scattered: Ray,
    pub albedo: Color,
}

impl Scatter {
    pub fn new(scattered: Ray, albedo: Color) -> Scatter {
        Scatter { scattered, albedo }
    }
}

pub trait Material {
    fn scatter(&self, ray: Ray, rec: &HitRecord, rng: &mut ThreadRng) -> Option<Scatter>;
}

pub struct Lambertian<T: Texture> {
    pub albedo: T,
}

impl<T: Texture> Lambertian<T> {
    pub fn new(albedo: T) -> Lambertian<T> {
        Lambertian { albedo }
    }
}

impl<T: Texture> Material for Lambertian<T> {
    fn scatter(&self, ray: Ray, rec: &HitRecord, rng: &mut ThreadRng) -> Option<Scatter> {
        let scatter_direction = rec.normal + Vec3::random_unit_vector(rng);
        let scattered = Ray::new(rec.point, scatter_direction, ray.time);
        Some(Scatter::new(
            scattered,
            self.albedo.value(rec.u, rec.v, &rec.point),
        ))
    }
}

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
    fn scatter(&self, ray: Ray, rec: &HitRecord, rng: &mut ThreadRng) -> Option<Scatter> {
        let reflected = ray.direction.unit().reflect(rec.normal);
        let scattered = Ray::new(
            rec.point,
            reflected + self.fuzz * Vec3::random_in_unit_sphere(rng),
            ray.time,
        );
        if scattered.direction.dot(&rec.normal) > 0.0 {
            Some(Scatter::new(scattered, self.albedo))
        } else {
            None
        }
    }
}

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
    fn scatter(&self, ray: Ray, rec: &HitRecord, rng: &mut ThreadRng) -> Option<Scatter> {
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
            Scatter::new(Ray::new(rec.point, reflected, ray.time), attenuation)
        } else if utils::random_double(rng) < Dielectric::schlick(cos_theta, etai_over_etat) {
            let reflected = unit_direction.reflect(rec.normal);
            Scatter::new(Ray::new(rec.point, reflected, ray.time), attenuation)
        } else {
            let refracted = unit_direction.refract(rec.normal, etai_over_etat);
            Scatter::new(Ray::new(rec.point, refracted, ray.time), attenuation)
        };
        Some(scattered)
    }
}
