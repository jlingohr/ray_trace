use ray_trace::ray::Ray;
use ray_trace::vector::Vec3;
use ray_trace::color::Color;
use ray_trace::point::Point;
use ray_trace::hittable;
use ray_trace::sphere::Sphere;
use std::rc::Rc;


fn main() {
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: u32 = 384;
    const IMAGE_HEIGHT: u32 = ((IMAGE_WIDTH as f64) / ASPECT_RATIO) as u32;

    println!("P3\n{} {}\n255", IMAGE_WIDTH, IMAGE_HEIGHT);

    let viewport_height = 2.0;
    let viewport_width = ASPECT_RATIO * viewport_height;
    let focal_length = 1.0;

    let origin = Point::new(0.0, 0.0, 0.0);
    let horizontal = Vec3::new(viewport_width, 0.0, 0.0);
    let vertical = Vec3::new(0.0, viewport_height, 0.0);
    let lower_left_corner = origin - (horizontal/2.0) - (vertical/2.0) - Vec3::new(0.0, 0.0, focal_length);


    let mut world = hittable::HittableList::new();
    world.add(Rc::new(Sphere::new(Point::new(0.0, 0.0, -1.0), 0.5)));
    world.add(Rc::new(Sphere::new(Point::new(0.0, -100.5, -1.0), 100.0)));

    for j in (0..IMAGE_HEIGHT).rev(){
        for i in 0..IMAGE_WIDTH {
            let u = (i as f64) / (IMAGE_WIDTH-1) as f64;
            let v = (j as f64) / (IMAGE_HEIGHT-1) as f64;
            let r = Ray::new(origin, lower_left_corner + u*horizontal + v*vertical - origin);
            let pixel_color = ray_color(r, &world);
            pixel_color.write_color();
        }
    }
}

fn ray_color(r: Ray, world: &dyn hittable::Hittable) -> Color {
    if let Some(hit) = world.hit(&r, 0.0, f64::INFINITY) {
        return 0.5 * (hit.normal + Color::new(1.0, 1.0, 1.0))
    }
    let unit_direction = r.direction.unit();
    let t = 0.5 * (unit_direction.y + 1.0);
    return ((1.0 - t) * Color::new(1.0, 1.0, 1.0)) + t * Color::new(0.5, 0.7, 1.0)
}


