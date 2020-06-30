use ray_trace::camera::Camera;
use ray_trace::point::Point;
use ray_trace::run;
use ray_trace::scenes;
use ray_trace::vector::Vec3;
use ray_trace::Config;
use ray_trace::Image;

use std::env;
use std::process;

fn main() {
    const ASPECT_RATIO: f64 = 1.0 / 1.0;
    const IMAGE_WIDTH: u32 = 600;
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
    // let world = scenes::two_perlin_spheres(&mut rng);
    // let world = scenes::earth(&mut rng);
    // let world = scenes::simple_light(&mut rng);
    // let world = scenes::cornell_box(&mut rng);
    // let world = scenes::cornell_smoke(&mut rng);
    let world = scenes::final_scene(&mut rng);

    let look_from = Point::new(478.0, 278.0, -600.0);
    let look_at = Point::new(278.0, 278.0, 0.0);
    let view_up = Vec3::new(0.0, 1.0, 0.0);
    let dist_to_focus = 10.0;
    let aperature = 0.0;
    let vfov = 40.0;
    let cam = Camera::new(
        look_from,
        look_at,
        view_up,
        vfov,
        ASPECT_RATIO,
        aperature,
        dist_to_focus,
        0.0,
        1.0,
    );

    let image = run(
        cam,
        world,
        IMAGE_HEIGHT,
        IMAGE_WIDTH,
        SAMPLES_PER_PIXEL,
        MAX_DEPTH,
    );
    Image::write(
        config.filename,
        image,
        IMAGE_WIDTH,
        IMAGE_HEIGHT,
        SAMPLES_PER_PIXEL,
    );
}
