use ray_trace::camera::Camera;
use ray_trace::color::Color;
use ray_trace::hittable;
use ray_trace::point::Point;
use ray_trace::ray::Ray;
use ray_trace::scenes;
use ray_trace::utils;
use ray_trace::vector::Vec3;

use rand::prelude::ThreadRng;
use std::env;
use std::fs;
use std::process;

fn main() {
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: u32 = 384;
    const IMAGE_HEIGHT: u32 = ((IMAGE_WIDTH as f64) / ASPECT_RATIO) as u32;
    const SAMPLES_PER_PIXEL: u32 = 100;
    const MAX_DEPTH: u32 = 50;

    let args: Vec<String> = env::args().collect();
    let config = Config::new(&args).unwrap_or_else(|err| {
        println!("Problem parsing arguments: {}", err);
        process::exit(1);
    });

    let mut rng = rand::thread_rng();
    // let world = scenes::random_scene(&mut rng);
    // let world = scenes::random_bvh_scene(&mut rng);
    // let world = scenes::random_checkered_scene(&mut rng);
    // let world = scenes::two_spheres(&mut rng);
    let world = scenes::two_perlin_spheres(&mut rng);

    let look_from = Point::new(13.0, 2.0, 3.0);
    let look_at = Point::new(0.0, 0.0, 0.0);
    let view_up = Vec3::new(0.0, 1.0, 0.0);
    let dist_to_focus = 10.0;
    let aperature = 0.0;
    let cam = Camera::new(
        look_from,
        look_at,
        view_up,
        20.0,
        ASPECT_RATIO,
        aperature,
        dist_to_focus,
        0.0,
        1.0,
    );

    let mut image_str: Vec<String> = Vec::new();
    image_str.push(format!("P3\n{} {}\n255\n", IMAGE_WIDTH, IMAGE_HEIGHT));

    for j in (0..IMAGE_HEIGHT).rev() {
        println!("\rScanlines remaining: {}", j);
        for i in 0..IMAGE_WIDTH {
            let mut pixel_color = Color::new(0.0, 0.0, 0.0);
            for _s in 0..SAMPLES_PER_PIXEL {
                let u = ((i as f64) + utils::random_double(&mut rng)) / (IMAGE_WIDTH - 1) as f64;
                let v = ((j as f64) + utils::random_double(&mut rng)) / (IMAGE_HEIGHT - 1) as f64;
                let r = cam.get_ray(&mut rng, u, v);
                pixel_color += ray_color(r, &world, MAX_DEPTH, &mut rng);
            }
            image_str.push(pixel_color.write_color(SAMPLES_PER_PIXEL));
        }
    }
    let image = image_str.join("");
    if fs::write(config.filename, image).is_err() {
        eprintln!("Could not generate image");
    }
}

struct Config {
    filename: String,
}

impl Config {
    fn new(args: &[String]) -> Result<Config, &'static str> {
        if args.len() < 2 {
            return Err("Not enough arguments");
        }
        let filename = args[1].clone();
        Ok(Config { filename })
    }
}

fn ray_color(
    r: Ray,
    world: &Box<dyn hittable::Hittable>,
    depth: u32,
    rng: &mut ThreadRng,
) -> Color {
    if depth <= 0 {
        return Color::new(0.0, 0.0, 0.0);
    }
    if let Some(hit) = world.hit(&r, 0.001, f64::INFINITY) {
        if let Some(scattered) = hit.material.scatter(r, &hit, rng) {
            return scattered.albedo * ray_color(scattered.scattered, world, depth - 1, rng);
        }
        return Color::new(0.0, 0.0, 0.0);
    }
    let unit_direction = r.direction.unit();
    let t = 0.5 * (unit_direction.y + 1.0);
    return ((1.0 - t) * Color::new(1.0, 1.0, 1.0)) + t * Color::new(0.5, 0.7, 1.0);
}
