use ray_trace::ray::Ray;
use ray_trace::vector::Vec3;
use ray_trace::color::Color;
use ray_trace::point::Point;
use ray_trace::hittable;
use ray_trace::sphere::Sphere;
use ray_trace::camera::Camera;
use ray_trace::utils;
use ray_trace::material::{Lambertian, Metal, Dielectric};

use std::rc::Rc;
use rand::Rng;
use rand::prelude::ThreadRng;

fn main() {
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: u32 = 384;
    const IMAGE_HEIGHT: u32 = ((IMAGE_WIDTH as f64) / ASPECT_RATIO) as u32;
    const SAMPLES_PER_PIXEL: u32 = 100;
    const MAX_DEPTH: u32 = 50;

    println!("P3\n{} {}\n255", IMAGE_WIDTH, IMAGE_HEIGHT);

    let mut world = hittable::HittableList::new();
    world.add(Rc::new(Sphere::new(Point::new(0.0, 0.0, -1.0), 0.5, Rc::new(Lambertian::new(Color::new(0.1, 0.2, 0.5))))));
    world.add(Rc::new(Sphere::new(Point::new(0.0, -100.5, -1.0), 100.0, Rc::new(Lambertian::new(Color::new(0.8, 0.8, 0.0))))));
    world.add(Rc::new(Sphere::new(Point::new(1.0, 0.0, -1.0), 0.5, Rc::new(Metal::new(Color::new(0.8, 0.6, 0.2), 0.0)))));
    world.add(Rc::new(Sphere::new(Point::new(-1.0, 0.0, -1.0), 0.5, Rc::new(Dielectric::new(1.5)))));

    let cam = Camera::new();
    let mut rng = rand::thread_rng();

    for j in (0..IMAGE_HEIGHT).rev(){
        for i in 0..IMAGE_WIDTH {
            let mut pixel_color = Color::new(0.0, 0.0, 0.0);
            for _s in 0..SAMPLES_PER_PIXEL {
                let u = ((i as f64) + utils::random_double(&mut rng)) / (IMAGE_WIDTH-1) as f64;
                let v = ((j as f64) + utils::random_double(&mut rng)) / (IMAGE_HEIGHT-1) as f64;
                let r = cam.get_ray(u, v);
                pixel_color += ray_color(r, &world, MAX_DEPTH, &mut rng);
            }
            pixel_color.write_color(SAMPLES_PER_PIXEL);
        }
    }
}

fn ray_color(r: Ray, world: &dyn hittable::Hittable, depth: u32, rng: &mut ThreadRng) -> Color {
    if depth <= 0 {
        return Color::new(0.0, 0.0, 0.0)
    }
    if let Some(hit) = world.hit(&r, 0.001, f64::INFINITY) {
        // let attenuation = Color::new(1.0, 1.0, 1.0);
        if let Some(scattered) = hit.material.scatter(r, &hit, rng) {
            return scattered.albedo * ray_color(scattered.scattered, world, depth-1, rng)
        }
        return Color::new(0.0, 0.0, 0.0)
    }
    let unit_direction = r.direction.unit();
    let t = 0.5 * (unit_direction.y + 1.0);
    return ((1.0 - t) * Color::new(1.0, 1.0, 1.0)) + t * Color::new(0.5, 0.7, 1.0)
}


