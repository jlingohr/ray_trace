use super::hittable::Hittable;
use super::onb::ONB;
use super::point::Point;
use super::vector::Vec3;
use rand::Rng;

pub fn random_cosine_direction() -> Vec3 {
    let mut rng = rand::thread_rng();
    let r1 = rng.gen::<f64>();
    let r2 = rng.gen::<f64>();
    let z = (1.0 - r2).sqrt();

    let phi = 2.0 * std::f64::consts::PI * r1;
    let x = phi.cos() * r2.sqrt();
    let y = phi.sin() * r2.sqrt();

    Vec3::new(x, y, z)
}

pub fn random_to_sphere(radius: f64, distance_squared: f64) -> Vec3 {
    let mut rng = rand::thread_rng();
    let r1 = rng.gen::<f64>();
    let r2 = rng.gen::<f64>();
    let z = 1.0 + r2 * (1.0 - radius * radius / distance_squared).sqrt() - 1.0;

    let phi = 2.0 * std::f64::consts::PI * r1;
    let x = phi.cos() * (1.0 - z * z).sqrt();
    let y = phi.sin() * (1.0 - z - z).sqrt();

    Vec3::new(x, y, z)
}

pub enum PDF<'a> {
    Cosine {
        uvw: ONB,
    },
    Hittable {
        origin: Point,
        hittable: &'a Box<dyn Hittable>,
    },
    Mixture {
        p: &'a PDF<'a>,
        q: &'a PDF<'a>,
    },
}

impl<'a> PDF<'a> {
    pub fn cosine_pdf(w: &Vec3) -> PDF<'a> {
        PDF::Cosine {
            uvw: ONB::build_from_w(*w),
        }
    }

    pub fn hittable_pdf(hittable: &'a Box<dyn Hittable>, origin: Point) -> PDF<'a> {
        PDF::Hittable { origin, hittable }
    }

    pub fn mixture_pdf(p: &'a PDF, q: &'a PDF) -> PDF<'a> {
        PDF::Mixture { p, q }
    }

    pub fn value(&self, direction: Vec3) -> f64 {
        match self {
            PDF::Cosine { uvw } => {
                let cosine = direction.unit().dot(&uvw.w());
                if cosine <= 0.0 {
                    1.0
                } else {
                    cosine / std::f64::consts::PI
                }
            }
            PDF::Hittable { origin, hittable } => hittable.pdf_value(*origin, direction),
            PDF::Mixture { p, q } => 0.5 * p.value(direction) + 0.5 * q.value(direction),
        }
    }

    pub fn generate(&self) -> Vec3 {
        match self {
            PDF::Cosine { uvw } => uvw.local(random_cosine_direction()),
            PDF::Hittable { origin, hittable } => hittable.random(*origin),
            PDF::Mixture { p, q } => {
                let mut rng = rand::thread_rng();
                if rng.gen::<bool>() {
                    p.generate()
                } else {
                    q.generate()
                }
            }
        }
    }
}
