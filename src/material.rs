use super::ray::Ray;
use super::hittable::HitRecord;
use super::color::Color;
use super::vector::Vec3;
use rand::prelude::ThreadRng;

pub struct Scatter {
    pub scattered: Ray,
    pub albedo: Color
}

impl Scatter {
    pub fn new(scattered: Ray, albedo: Color) -> Scatter {
        Scatter { scattered, albedo, }
    }
}

pub trait Material {
    fn scatter(&self, ray: Ray, rec: &HitRecord, rng: &mut ThreadRng) -> Option<Scatter>;
}

pub struct Lambertian {
    pub albedo: Color,
}

impl Lambertian {
    pub fn new(albedo: Color) -> Lambertian {
        Lambertian { albedo }
    }
}

impl Material for Lambertian {
    fn scatter(&self, _: Ray, rec: &HitRecord, rng: &mut ThreadRng) -> Option<Scatter> {
        let scatter_direction = rec.normal + Vec3::random_unit_vector(rng);
        let scattered = Ray::new(rec.point, scatter_direction);
        Some(Scatter::new(scattered, self.albedo))
    }
}

pub struct Metal {
    pub albedo: Color,
}

impl Metal {
    pub fn new(albedo: Color) -> Metal {
        Metal { albedo }
    }
}

impl Material for Metal {
    fn scatter(&self, ray: Ray, rec: &HitRecord, _rng: &mut ThreadRng) -> Option<Scatter> {
        let reflected = ray.direction.unit().reflect(rec.normal);
        let scattered = Ray::new(rec.point, reflected);
        if scattered.direction.dot(&rec.normal) > 0.0 {
            Some(Scatter::new(scattered, self.albedo))
        } else {
            None
        }
    }
}