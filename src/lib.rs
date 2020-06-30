pub mod aabb;
pub mod aarect;
pub mod bvh;
pub mod camera;
pub mod color;
pub mod cube;
pub mod hittable;
pub mod material;
pub mod medium;
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
use crate::ray::Ray;

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
    mut r: Ray,
    world: &Box<dyn hittable::Hittable>,
    depth: u32,
    rng: &mut ThreadRng,
) -> Color {
    let mut acc = Color::new(0.0, 0.0, 0.0);
    let mut strength = Color::new(1.0, 1.0, 1.0);
    let mut bounces = 0;

    while let Some(hit) = world.hit(&r, 0.001, std::f64::INFINITY) {
        acc += strength * hit.material.emitted(hit.u, hit.t, &hit.point);

        if let Some(scattered) = hit.material.scatter(r, &hit, rng) {
            r = scattered.scattered;
            strength = strength * scattered.albedo;
        } else {
            return acc;
        }

        if bounces == depth {
            return acc;
        }

        bounces += 1;
    }
    Color::new(0.0, 0.0, 0.0)
}

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
                ray_color(r, &world, max_depth, &mut rng)
            })
            .sum();
        col
    })
}
