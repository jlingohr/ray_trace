use std::cmp::min;
use std::prelude::v1::Iterator;

use nalgebra::Point2;

// use crate::camera::Camera;
use crate::camera::Camera;
use crate::geometry::bounds::Bounds2i;
use crate::geometry::ray::Ray;
use crate::integrators::whitted::WhittedIntegrator;
use crate::pbrt::{Float, Spectrum};
use crate::samplers::uniform::UniformSampler;
use crate::scene::Scene;

struct Point2Iterator {
    point: Point2<i32>,
    x_max: i32,
    y_max: i32,
}

impl Iterator for Point2Iterator {
    type Item = Point2<i32>;

    fn next(&mut self) -> Option<Self::Item> {
        self.point.x += 1;
        if self.point.x == self.x_max {
            self.point.x = 0;
            self.point.y += 1;
        }

        if self.point.y == self.y_max {
            None
        } else {
            Some(Point2::new(self.point.x, self.point.y))
        }
    }
}

// impl IntoIterator for Point2<i32> {
//     type Item = Point2<i32>;
//     type IntoIter = Point2Iterator;
//
//     fn into_iter(self) -> Self::IntoIter {
//         Point2Iterator {
//             point: Point2::new(-1, 0),
//             x_max: self.x,
//             y_max: self.y,
//         }
//     }
// }

pub enum SamplerIntegrator {
    Whitted(WhittedIntegrator),
}

//TODO make generic sampler
impl SamplerIntegrator {
    fn preprocess(&self, scene: &Scene, sampler: &mut UniformSampler) {
        unimplemented!()
    }

    pub fn render(&self, scene: &Scene) {
        let camera = self.get_camera();
        let film = camera.get_film();
        let sample_bounds = film.get_sample_bounds();
        let sample_extent = sample_bounds.diagonal();
        let tile_size = 16; //TODO don't hardcode this
        let n_tiles = Point2::new(
            (sample_extent.x + tile_size - 1) / tile_size,
            (sample_extent.y + tile_size - 1) / tile_size,
        );

        let sampler = self.get_sampler();

        let iterator = Point2Iterator {
            point: Point2::new(-1, 0),
            x_max: n_tiles.x,
            y_max: n_tiles.y,
        };

        for tile in iterator {
            //TODO parallel for
            // allocate memory for tile
            // let arena: MemoryArena;

            // get sampler instance for tile
            let seed = tile.y * n_tiles.x + tile.x;
            let mut tile_sampler = sampler.clone_with_seed(seed as u64);

            // compute sample bounds for tile
            let x0 = (sample_bounds.p_min.x + tile.x * tile_size) as i32;
            let x1 = min(x0 + tile_size as i32, sample_bounds.p_max.x as i32);
            let y0 = (sample_bounds.p_min.y + tile.y * tile_size) as i32;
            let y1 = min(y0 + tile_size, sample_bounds.p_max.y as i32);
            let tile_bounds = Bounds2i::new(Point2::new(x0, y0), Point2::new(x1, y1));

            // get film tile for tile
            let mut film_tile = film.get_film_tile(&tile_bounds);

            //loop over pixels in tile to render them
            for pixel in tile_bounds {
                // tile_sampler.start_pixel(&pixel);
                for _ in 0..tile_sampler.samples_per_pixel {
                    //TODO might need interior mutability
                    // initialize camera sample for current sample
                    // let camera_sample = tile_sampler.get_camera_sample(&pixel);
                    let camera_sample = tile_sampler.sample(&pixel);

                    // generate camera ray for current sample
                    let (ray, ray_weight) = camera.generate_ray_differential(&camera_sample);
                    let scaled_ray = ray
                        .scale_differential(1.0 / (tile_sampler.samples_per_pixel as Float).sqrt());

                    // evaluate radiance along camera ray
                    let lighting = if ray_weight >= 0.0 {
                        //TODO warnings and invalid radiance?
                        self.li(&scaled_ray, scene, &mut tile_sampler)
                    // issue warning if unexpected radiance value is returned
                    } else {
                        Spectrum::new(0.0)
                    };

                    // add camera rays contribution to image
                    film_tile.add_sample(&camera_sample.p_film, &lighting, ray_weight);
                }
            }
            // merge image into film
            film.merge_film_tile(&film_tile);
            // camera.get_film().merge_film_tile(&film_tile);
        }

        // Save final image after rendering
        camera.save_image();
    }

    fn get_camera(&self) -> &Camera {
        match self {
            SamplerIntegrator::Whitted(integrator) => integrator.get_camera(),
        }
    }

    fn get_sampler(&self) -> &UniformSampler {
        match self {
            SamplerIntegrator::Whitted(integrator) => integrator.get_sampler(),
        }
    }

    fn li(&self, ray: &Ray, scene: &Scene, sampler: &mut UniformSampler) -> Spectrum {
        match self {
            SamplerIntegrator::Whitted(integrator) => integrator.li(ray, scene, sampler, 0),
        }
    }
}
