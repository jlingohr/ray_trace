pub mod aabb;
pub mod aarect;
pub mod bvh;
pub mod camera;
pub mod color;
pub mod cube;
pub mod hittable;
pub mod material;
pub mod medium;
pub mod onb;
pub mod pdf;
pub mod perlin;
pub mod point;
pub mod ray;
pub mod rotate;
pub mod scenes;
pub mod sphere;
pub mod texture;
pub mod utils;
pub mod vector;

use crate::camera::Camera;
use crate::color::Color;
use crate::hittable::Hittable;
use crate::material::ScatterRecord;
use crate::pdf::PDF;
use crate::ray::Ray;
use crate::scenes::Scene;

use rand::prelude::ThreadRng;
use rayon::prelude::*;

use std::fs;

pub struct Config {
    pub filename: String,
}

impl Config {
    pub fn new(args: &[String]) -> Result<Config, &'static str> {
        if args.len() < 2 {
            return Err("Not enough arguments");
        }
        let filename = args[1].clone();
        Ok(Config { filename })
    }
}

pub fn ray_color(
    ray_in: Ray,
    world: &Box<dyn hittable::Hittable>,
    lights: &Box<dyn hittable::Hittable>,
    depth: u32,
    rng: &mut ThreadRng,
) -> Color {
    if depth <= 0 {
        return Color::new(0.0, 0.0, 0.0);
    }

    if let Some(hit) = world.hit(&ray_in, 0.001, std::f64::INFINITY) {
        let emitted = hit.material.emitted(ray_in, &hit);

        if let Some(scatter_record) = hit.material.scatter(ray_in, &hit, rng) {
            match scatter_record {
                ScatterRecord::Specular {
                    scattered_ray,
                    albedo,
                } => return albedo * ray_color(scattered_ray, world, lights, depth - 1, rng),
                ScatterRecord::Scatter { albedo, pdf } => {
                    let lights_ptr = PDF::hittable_pdf(&lights, hit.point);
                    let p = PDF::mixture_pdf(&lights_ptr, &pdf);
                    let scattered = Ray::new(hit.point, p.generate(), ray_in.time);
                    let pdf_val = p.value(scattered.direction);
                    return emitted
                        + albedo
                            * hit.material.scattering_pdf(&ray_in, &hit, &scattered)
                            * ray_color(scattered, world, lights, depth - 1, rng)
                            / pdf_val;
                }
            }
        } else {
            return emitted;
        }
    } else {
        return Color::new(0.0, 0.0, 0.0);
    }
}

// pub fn ray_color(
//     mut r: Ray,
//     world: &Box<dyn hittable::Hittable>,
//     depth: u32,
//     rng: &mut ThreadRng,
// ) -> Color {
//     let mut acc = Color::new(0.0, 0.0, 0.0);
//     let mut strength = Color::new(1.0, 1.0, 1.0);
//     let mut bounces = 0;

//     while let Some(hit) = world.hit(&r, 0.001, std::f64::INFINITY) {
//         acc += strength * hit.material.emitted(&hit);

//         if let Some(scatter_record) = hit.material.scatter(r, &hit, rng) {
//             match scatter_record {
//                 ScatterRecord::Specular { ray, albedo } => {
//                     r = ray;
//                     strength = strength * albedo;
//                 }
//                 ScatterRecord::Scatter { ray, pdf, albedo } => {
//                     r = ray;
//                     strength = strength * albedo;
//                 }
//             }
//         } else {
//             return acc;
//         }

//         if bounces == depth {
//             return acc;
//         }

//         bounces += 1;
//     }
//     Color::new(0.0, 0.0, 0.0)
// }

pub struct Image(Vec<Vec<Color>>);

impl Image {
    pub fn compute(nx: u32, ny: u32, f: impl Fn(u32, u32) -> Color + Sync) -> Image {
        Image(
            (0..ny)
                .into_par_iter()
                .rev()
                .map(|y| (0..nx).map(|x| f(x, y)).collect())
                .collect(),
        )
    }

    pub fn write(filename: String, image: Image, width: u32, height: u32, num_samples: u32) {
        let mut image_str: Vec<String> = Vec::new();
        image_str.push(format!("P3\n{} {}\n255\n", width, height));
        for scanline in image.0 {
            for col in scanline {
                image_str.push(col.write_color(num_samples));
            }
        }

        let image_str = image_str.join("");
        if fs::write(filename, image_str).is_err() {
            eprintln!("Could not generate image");
        }
    }
}

pub fn run(
    cam: Camera,
    world: Box<dyn Hittable>,
    lights: Box<dyn Hittable>,
    height: u32,
    width: u32,
    num_samples: u32,
    max_depth: u32,
) -> Image {
    Image::compute(height, width, |x, y| {
        let col: Color = (0..num_samples)
            .map(|_| {
                let mut rng = rand::thread_rng();
                let u = ((x as f64) + utils::random_double(&mut rng)) / (width - 1) as f64;
                let v = ((y as f64) + utils::random_double(&mut rng)) / (height - 1) as f64;
                let r = cam.get_ray(&mut rng, u, v);
                ray_color(r, &world, &lights, max_depth, &mut rng)
            })
            .sum();
        col
    })
}
