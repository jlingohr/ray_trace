use super::aarect::{AARectangle, Plane};
use super::bvh::BVH;
use super::color::Color;
use super::hittable::{FlipFace, Hittable, HittableList};
use super::material::{Dielectric, DiffuseLight, Lambertian, Metal};
use super::point::Point;
use super::sphere::{MovingSphere, Sphere};
use super::texture::{Checkered, ImageTexture, NoiseTexture, SolidColor};
use super::utils;
use super::vector::Vec3;

use image;
use rand::prelude::ThreadRng;

pub fn random_scene(rng: &mut ThreadRng) -> Box<dyn Hittable> {
    let mut world = HittableList::new();
    let ground_material = Lambertian::new(SolidColor::new(0.5, 0.5, 0.5));
    world.add(Sphere::new(
        Point::new(0.0, -1000.0, 0.0),
        1000.0,
        ground_material,
    ));

    for a in -11..11 {
        for b in -11..11 {
            let choose_mat = utils::random_double(rng);
            let center = Point::new(
                (a as f64) + (0.9 * utils::random_double(rng)),
                0.2,
                (b as f64) + (0.9 * utils::random_double(rng)),
            );
            if (center - Point::new(4.0, 0.2, 0.0)).len() > 0.9 {
                if choose_mat < 0.8 {
                    // diffuse
                    let albedo = Color::random(rng) * Color::random(rng);
                    let material = Lambertian::new(SolidColor::new(albedo.r, albedo.g, albedo.b));
                    let center2 =
                        center + Vec3::new(0.0, utils::random_in_range(rng, 0.0, 0.5), 0.0);
                    world.add(MovingSphere::new(center, center2, 0.0, 1.0, 0.2, material));
                } else if choose_mat < 0.95 {
                    // metal
                    let albedo = Color::random_in_range(rng, 0.5, 1.0);
                    let fuzz = utils::random_in_range(rng, 0.0, 0.5);
                    let material = Metal::new(albedo, fuzz);
                    world.add(Sphere::new(center, 0.2, material));
                } else {
                    // glass
                    let material = Dielectric::new(1.5);
                    world.add(Sphere::new(center, 0.2, material));
                };
            }
        }
    }

    let material1 = Dielectric::new(1.5);
    world.add(Sphere::new(Point::new(0.0, 1.0, 0.0), 1.0, material1));

    let material2 = Lambertian::new(SolidColor::new(0.4, 0.2, 0.1));
    world.add(Sphere::new(Point::new(-4.0, 1.0, 0.0), 1.0, material2));

    let material3 = Metal::new(Color::new(0.7, 0.6, 0.5), 0.0);
    world.add(Sphere::new(Point::new(4.0, 1.0, 0.0), 1.0, material3));

    Box::new(world)
}

pub fn random_checkered_scene(rng: &mut ThreadRng) -> Box<dyn Hittable> {
    let mut world = HittableList::new();
    let checkered = Checkered::new(
        SolidColor::new(0.2, 0.3, 0.1),
        SolidColor::new(0.9, 0.9, 0.9),
    );
    let ground_material = Lambertian::new(checkered);
    world.add(Sphere::new(
        Point::new(0.0, -1000.0, 0.0),
        1000.0,
        ground_material,
    ));

    for a in -11..11 {
        for b in -11..11 {
            let choose_mat = utils::random_double(rng);
            let center = Point::new(
                (a as f64) + (0.9 * utils::random_double(rng)),
                0.2,
                (b as f64) + (0.9 * utils::random_double(rng)),
            );
            if (center - Point::new(4.0, 0.2, 0.0)).len() > 0.9 {
                if choose_mat < 0.8 {
                    // diffuse
                    let albedo = Color::random(rng) * Color::random(rng);
                    let material = Lambertian::new(SolidColor::new(albedo.r, albedo.g, albedo.b));
                    let center2 =
                        center + Vec3::new(0.0, utils::random_in_range(rng, 0.0, 0.5), 0.0);
                    world.add(MovingSphere::new(center, center2, 0.0, 1.0, 0.2, material));
                } else if choose_mat < 0.95 {
                    // metal
                    let albedo = Color::random_in_range(rng, 0.5, 1.0);
                    let fuzz = utils::random_in_range(rng, 0.0, 0.5);
                    let material = Metal::new(albedo, fuzz);
                    world.add(Sphere::new(center, 0.2, material));
                } else {
                    // glass
                    let material = Dielectric::new(1.5);
                    world.add(Sphere::new(center, 0.2, material));
                };
            }
        }
    }

    let material1 = Dielectric::new(1.5);
    world.add(Sphere::new(Point::new(0.0, 1.0, 0.0), 1.0, material1));

    let material2 = Lambertian::new(SolidColor::new(0.4, 0.2, 0.1));
    world.add(Sphere::new(Point::new(-4.0, 1.0, 0.0), 1.0, material2));

    let material3 = Metal::new(Color::new(0.7, 0.6, 0.5), 0.0);
    world.add(Sphere::new(Point::new(4.0, 1.0, 0.0), 1.0, material3));

    Box::new(world)
}

pub fn two_spheres(rng: &mut ThreadRng) -> Box<dyn Hittable> {
    let mut world = HittableList::new();
    let checkered = Checkered::new(
        SolidColor::new(0.2, 0.3, 0.1),
        SolidColor::new(0.9, 0.9, 0.9),
    );
    world.add(Sphere::new(
        Point::new(0.0, -10.0, 0.0),
        10.0,
        Lambertian::new(checkered.clone()),
    ));
    world.add(Sphere::new(
        Point::new(0.0, 10.0, 0.0),
        10.0,
        Lambertian::new(checkered),
    ));

    Box::new(world)
}

pub fn two_perlin_spheres(rng: &mut ThreadRng) -> Box<dyn Hittable> {
    let mut world = HittableList::new();
    let pertext = NoiseTexture::new(1.0);
    world.add(Sphere::new(
        Point::new(0.0, -1000.0, 0.0),
        1000.0,
        Lambertian::new(pertext.clone()),
    ));
    world.add(Sphere::new(
        Point::new(0.0, 2.0, 0.0),
        2.0,
        Lambertian::new(pertext),
    ));

    Box::new(world)
}

pub fn earth(rng: &mut ThreadRng) -> Box<dyn Hittable> {
    let img = image::open("earthmap.jpg")
        .expect("image not found")
        .to_rgb();
    let (nx, ny) = img.dimensions();
    let data = img.into_raw();
    let mut world = HittableList::new();
    let earth_texture = ImageTexture::new(data, nx, ny);
    let earth_surface = Lambertian::new(earth_texture);
    world.add(Sphere::new(Point::new(0.0, 0.0, 0.0), 2.0, earth_surface));

    Box::new(world)
}

pub fn simple_light(rng: &mut ThreadRng) -> Box<dyn Hittable> {
    let mut world = HittableList::new();
    let pertext = NoiseTexture::new(4.0);
    world.add(Sphere::new(
        Point::new(0.0, -1000.0, 0.0),
        1000.0,
        Lambertian::new(pertext.clone()),
    ));
    world.add(Sphere::new(
        Point::new(0.0, 2.0, 0.0),
        2.0,
        Lambertian::new(pertext),
    ));

    let diff_light = DiffuseLight::new(SolidColor::new(4.0, 4.0, 4.0));
    world.add(Sphere::new(
        Point::new(0.0, 7.0, 0.0),
        2.0,
        diff_light.clone(),
    ));
    world.add(AARectangle::new(
        3.0,
        5.0,
        1.0,
        3.0,
        -2.0,
        diff_light,
        Plane::XY,
    ));
    Box::new(world)
}

pub fn cornell_box(rng: &mut ThreadRng) -> Box<dyn Hittable> {
    let mut world = HittableList::new();

    let red = Lambertian::new(SolidColor::new(0.65, 0.05, 0.05));
    let white = Lambertian::new(SolidColor::new(0.73, 0.73, 0.73));
    let green = Lambertian::new(SolidColor::new(0.12, 0.45, 0.15));
    let light = DiffuseLight::new(SolidColor::new(15.0, 15.0, 15.0));

    world.add(FlipFace::new(AARectangle::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        green,
        Plane::YZ,
    )));
    world.add(AARectangle::new(
        0.0,
        555.0,
        0.0,
        555.0,
        0.0,
        red,
        Plane::YZ,
    ));
    world.add(AARectangle::new(
        213.0,
        343.0,
        227.0,
        332.0,
        554.0,
        light,
        Plane::XZ,
    ));
    world.add(FlipFace::new(AARectangle::new(
        0.0,
        555.0,
        0.0,
        555.0,
        0.0,
        white.clone(),
        Plane::XZ,
    )));
    world.add(AARectangle::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        white.clone(),
        Plane::XZ,
    ));
    world.add(AARectangle::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        white,
        Plane::XY,
    ));

    Box::new(world)
}

pub fn random_bvh_scene(rng: &mut ThreadRng) -> Box<dyn Hittable> {
    let mut world = HittableList::new();
    let ground_material = Lambertian::new(SolidColor::new(0.5, 0.5, 0.5));
    world.add(Sphere::new(
        Point::new(0.0, -1000.0, 0.0),
        1000.0,
        ground_material,
    ));

    for a in -11..11 {
        for b in -11..11 {
            let choose_mat = utils::random_double(rng);
            let center = Point::new(
                (a as f64) + (0.9 * utils::random_double(rng)),
                0.2,
                (b as f64) + (0.9 * utils::random_double(rng)),
            );
            if (center - Point::new(4.0, 0.2, 0.0)).len() > 0.9 {
                if choose_mat < 0.8 {
                    // diffuse
                    let albedo = Color::random(rng) * Color::random(rng);
                    let material = Lambertian::new(SolidColor::new(albedo.r, albedo.g, albedo.b));
                    let center2 =
                        center + Vec3::new(0.0, utils::random_in_range(rng, 0.0, 0.5), 0.0);
                    world.add(MovingSphere::new(center, center2, 0.0, 1.0, 0.2, material));
                } else if choose_mat < 0.95 {
                    // metal
                    let albedo = Color::random_in_range(rng, 0.5, 1.0);
                    let fuzz = utils::random_in_range(rng, 0.0, 0.5);
                    let material = Metal::new(albedo, fuzz);
                    world.add(Sphere::new(center, 0.2, material));
                } else {
                    // glass
                    let material = Dielectric::new(1.5);
                    world.add(Sphere::new(center, 0.2, material));
                };
            }
        }
    }

    let material1 = Dielectric::new(1.5);
    world.add(Sphere::new(Point::new(0.0, 1.0, 0.0), 1.0, material1));

    let material2 = Lambertian::new(SolidColor::new(0.4, 0.2, 0.1));
    world.add(Sphere::new(Point::new(-4.0, 1.0, 0.0), 1.0, material2));

    let material3 = Metal::new(Color::new(0.7, 0.6, 0.5), 0.0);
    world.add(Sphere::new(Point::new(4.0, 1.0, 0.0), 1.0, material3));

    Box::new(BVH::new(world.objects, 0.0, 1.0))
}
