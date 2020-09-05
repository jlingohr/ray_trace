use std::fs;

use nalgebra::Point2;

use crate::geometry::bounds::Bounds2i;
use crate::pbrt::{clamp, Float};

// Write image to file
pub fn write_image(
    name: &String,
    rgb: Vec<Float>,
    output_bounds: &Bounds2i,
    total_resolution: &Point2<i32>,
) {
    let resolution = output_bounds.diagonal();
    write_image_ppm(
        name,
        rgb,
        resolution.x,
        resolution.y,
        total_resolution.x,
        total_resolution.y,
        output_bounds.p_min.x,
        output_bounds.p_min.y,
    )
}

// TODO ppm might be top-left, not bottom left
// TODO return Result<>
fn write_image_ppm(
    name: &String,
    pixels: Vec<Float>,
    x_res: i32,
    y_res: i32,
    width: i32,
    height: i32,
    x_offset: i32,
    y_offset: i32,
) {
    let mut image_str: Vec<String> = Vec::new();
    image_str.push(format!("P3\n{} {}\n255\n", width, height));

    for i in 0..x_res * y_res {
        let ir: u32 = (256.0 * clamp(pixels[3 * i as usize], 0.0, 0.999)) as u32;
        let ig: u32 = (256.0 * clamp(pixels[(3 * i + 1) as usize], 0.0, 0.999)) as u32;
        let ib: u32 = (256.0 * clamp(pixels[(3 * i + 2) as usize], 0.0, 0.999)) as u32;
        image_str.push(format!("{} {} {}\n", ir, ig, ib));
    }
    let image_str = image_str.join("");
    if fs::write(name, image_str).is_err() {
        eprintln!("Could not generate image");
    }
}
