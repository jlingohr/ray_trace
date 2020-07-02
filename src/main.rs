use ray_trace::run;
use ray_trace::scenes::Scene;
use ray_trace::Config;
use ray_trace::Image;

use std::env;
use std::process;

fn main() {
    const ASPECT_RATIO: f64 = 1.0 / 1.0;
    const IMAGE_WIDTH: u32 = 600;
    const IMAGE_HEIGHT: u32 = ((IMAGE_WIDTH as f64) / ASPECT_RATIO) as u32;
    const SAMPLES_PER_PIXEL: u32 = 1000;
    const MAX_DEPTH: u32 = 50;

    let args: Vec<String> = env::args().collect();
    let config = Config::new(&args).unwrap_or_else(|err| {
        println!("Problem parsing arguments: {}", err);
        process::exit(1);
    });

    let mut rng = rand::thread_rng();
    // let scene = Scene::final_scene(&mut rng, ASPECT_RATIO);
    let scene = Scene::cornell_box(&mut rng, ASPECT_RATIO);
    let cam = scene.cam;
    let world = scene.world;
    let lights = scene.lights;

    let image = run(
        cam,
        world,
        lights,
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
