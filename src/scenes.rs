use super::aarect::{AARectangle, Plane};
use super::bvh::BVH;
use super::color::Color;
use super::cube::Cube;
use super::hittable::{FlipFace, Hittable, HittableList, Translate};
use super::material::{Dielectric, DiffuseLight, Lambertian, Metal};
use super::medium::ConstantMedium;
use super::point::Point;
use super::rotate::RotateY;
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

pub fn two_spheres(_rng: &mut ThreadRng) -> Box<dyn Hittable> {
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

pub fn two_perlin_spheres(_rng: &mut ThreadRng) -> Box<dyn Hittable> {
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

pub fn earth(_rng: &mut ThreadRng) -> Box<dyn Hittable> {
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

pub fn simple_light(_rng: &mut ThreadRng) -> Box<dyn Hittable> {
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

pub fn cornell_box(_rng: &mut ThreadRng) -> Box<dyn Hittable> {
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
    let box1 = Cube::new(
        Point::new(0.0, 0.0, 0.0),
        Point::new(165.0, 330.0, 165.0),
        white.clone(),
    );
    let box1 = RotateY::new(box1, 15.0);
    let box1 = Translate::new(box1, Vec3::new(265.0, 0.0, 295.0));
    world.add(box1);

    let box2 = Cube::new(
        Point::new(0.0, 0.0, 0.0),
        Point::new(165.0, 165.0, 165.0),
        white.clone(),
    );
    let box2 = RotateY::new(box2, -18.0);
    let box2 = Translate::new(box2, Vec3::new(130.0, 0.0, 65.0));
    world.add(box2);

    Box::new(world)
}

pub fn cornell_smoke(_rng: &mut ThreadRng) -> Box<dyn Hittable> {
    let mut world = HittableList::new();

    let red = Lambertian::new(SolidColor::new(0.65, 0.05, 0.05));
    let white = Lambertian::new(SolidColor::new(0.73, 0.73, 0.73));
    let green = Lambertian::new(SolidColor::new(0.12, 0.45, 0.15));
    let light = DiffuseLight::new(SolidColor::new(7.0, 7.0, 7.0));

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
    let box1 = Cube::new(
        Point::new(0.0, 0.0, 0.0),
        Point::new(165.0, 330.0, 165.0),
        white.clone(),
    );
    let box1 = RotateY::new(box1, 15.0);
    let box1 = Translate::new(box1, Vec3::new(265.0, 0.0, 295.0));
    world.add(ConstantMedium::new(
        box1,
        0.01,
        SolidColor::new(0.0, 0.0, 0.0),
    ));

    let box2 = Cube::new(
        Point::new(0.0, 0.0, 0.0),
        Point::new(165.0, 165.0, 165.0),
        white.clone(),
    );
    let box2 = RotateY::new(box2, -18.0);
    let box2 = Translate::new(box2, Vec3::new(130.0, 0.0, 65.0));
    world.add(ConstantMedium::new(
        box2,
        0.01,
        SolidColor::new(1.0, 1.0, 1.0),
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

pub fn final_scene(rng: &mut ThreadRng) -> Box<dyn Hittable> {
    let mut boxes1 = HittableList::new();
    let ground = Lambertian::new(SolidColor::new(0.48, 0.83, 0.53));
    let boxes_per_side = 20;
    for i in 0..boxes_per_side {
        for j in 0..boxes_per_side {
            let w = 100.0;
            let x0 = -1000.0 + (i as f64) * w;
            let z0 = -1000.0 + (j as f64) * w;
            let y0 = 0.0;
            let x1 = x0 + w;
            let y1 = utils::random_in_range(rng, 1.0, 101.0);
            let z1 = z0 + w;

            boxes1.add(Cube::new(
                Point::new(x0, y0, z0),
                Point::new(x1, y1, z1),
                ground,
            ));
        }
    }

    let mut objects = HittableList::new();
    objects.add(BVH::new(boxes1.objects, 0.0, 1.0));
    let light = DiffuseLight::new(SolidColor::new(7.0, 7.0, 7.0));
    objects.add(AARectangle::new(
        123.0,
        423.0,
        147.0,
        412.0,
        554.0,
        light.clone(),
        Plane::XZ,
    ));

    let center1 = Point::new(400.0, 400.0, 200.0);
    let center2 = center1 + Vec3::new(30.0, 0.0, 0.0);
    let moving_sphere_material = Lambertian::new(SolidColor::new(0.7, 0.3, 0.1));
    objects.add(MovingSphere::new(
        center1,
        center2,
        0.0,
        1.0,
        50.0,
        moving_sphere_material,
    ));

    objects.add(Sphere::new(
        Point::new(260.0, 150.0, 45.0),
        70.0,
        Dielectric::new(1.5),
    ));
    objects.add(Sphere::new(
        Point::new(0.0, 150.0, 145.0),
        50.0,
        Metal::new(Color::new(0.8, 0.8, 0.9), 10.0),
    ));

    let boundary = Sphere::new(Point::new(360.0, 150.0, 145.0), 70.0, Dielectric::new(1.5));
    objects.add(boundary.clone());
    objects.add(ConstantMedium::new(
        boundary,
        0.2,
        SolidColor::new(0.2, 0.4, 0.9),
    ));

    let boundary = Sphere::new(Point::new(0.0, 0.0, 0.0), 5000.0, Dielectric::new(1.5));
    objects.add(ConstantMedium::new(
        boundary,
        0.0001,
        SolidColor::new(1.0, 1.0, 1.0),
    ));

    let img = image::open("earthmap.jpg")
        .expect("image not found")
        .to_rgb();
    let (nx, ny) = img.dimensions();
    let data = img.into_raw();
    let earth_texture = ImageTexture::new(data, nx, ny);
    let earth_surface = Lambertian::new(earth_texture);
    objects.add(Sphere::new(
        Point::new(400.0, 200.0, 400.0),
        100.0,
        earth_surface,
    ));

    let pertext = NoiseTexture::new(0.1);
    objects.add(Sphere::new(
        Point::new(220.0, 280.0, 300.0),
        80.0,
        Lambertian::new(pertext),
    ));

    let mut boxes2 = HittableList::new();
    let white = Lambertian::new(SolidColor::new(0.73, 0.73, 0.73));
    let ns = 1000;
    for j in 0..ns {
        boxes2.add(Sphere::new(
            Point::new(
                utils::random_in_range(rng, 0.0, 165.0),
                utils::random_in_range(rng, 0.0, 165.0),
                utils::random_in_range(rng, 0.0, 165.0),
            ),
            10.0,
            white,
        ));
    }

    objects.add(Translate::new(
        RotateY::new(BVH::new(boxes2.objects, 0.0, 1.0), 15.0),
        Vec3::new(-100.0, 270.0, 395.0),
    ));

    Box::new(objects)
}
