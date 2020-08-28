use crate::core::filter::Filter;
use crate::core::geometry::bounds::{Bounds2f, Bounds2i};
use crate::core::pbrt::{Float, Spectrum};
extern crate nalgebra as na;
use na::geometry::Point2;
use na::Vector2;
use std::cmp::{max, min};

const FILTER_TABLE_WIDTH: usize = 16;

// Represent individual pixels on film
// xyz: running weighted sums of spectral pixel contribution
// filter_weight_sum: sum of filter weight values for the sample contribution to the pixel
// splat_xyz: holds an unweighted sum of sample splats
// pad: unused. Purpose is to ensure Pixel structure is 32 or 64 byte size so minimize cache misses
pub struct Pixel {
    xyz: [Float; 3],
    filter_weight_sum: Float,
    splat_xyz: [Float; 3],
    pad: Float,
}

impl Pixel {
    pub fn new(
        lighting: Spectrum,
        filter_weight_sum: Float,
        splat_xyz: [Float; 3],
        pad: Float,
    ) -> Pixel {
        let xyz = lighting.to_xyz();
        Pixel {
            xyz,
            filter_weight_sum,
            splat_xyz,
            pad,
        }
    }
}

#[derive(Debug, Copy, Clone)]
struct FilmTilePixel {
    pub contrib_sum: Spectrum,
    pub filter_weight_sum: Float,
}

impl FilmTilePixel {
    pub fn new() -> FilmTilePixel {
        FilmTilePixel {
            contrib_sum: Spectrum::new(0.0),
            filter_weight_sum: 0.0,
        }
    }
}

#[derive(Debug, Copy, Clone)]
struct FilmTile<'a> {
    pixel_bounds: Bounds2i,
    filter_radius: Vector2<Float>,
    inv_filter_radius: Vector2<Float>,
    filter_table: &'a [Float; FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH],
    filter_table_size: usize,
    pixels: Vec<FilmTilePixel>,
}

impl<'a> FilmTile<'a> {
    pub fn new(
        pixel_bounds: Bounds2i,
        filter_radius: Vector2<Float>,
        inv_filter_radius: Vector2<Float>,
        filter_table: &'a [Float; FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH],
        filter_table_size: usize,
    ) -> FilmTile<'a> {
        let pixels: Vec<FilmTilePixel> =
            Vec::with_capacity(max(0, pixel_bounds.surface_area() as usize));

        FilmTile {
            pixel_bounds,
            filter_radius,
            inv_filter_radius,
            filter_table,
            filter_table_size,
            pixels,
        }
    }

    // Takes a sample and corresponding radiance value and weight and updates the stored image
    // using the reconstruction filter
    //TODO need to index using isize instead of i32
    pub fn add_sample(
        &mut self,
        p_film: &Point2<Float>,
        lighting: &Spectrum,
        sample_weight: Float,
    ) {
        // Compute sample's raster bounds
        let p_film_discrete = p_film - Vector2::new(0.5, 0.5);
        let mut p0: Point2<i32> = Point2::new(
            (p_film_discrete.x - &self.filter_radius).ceil(),
            (p_film_discrete.y - &self.filter_radius).ceil(),
        );
        let mut p1: Point2<i32> = Point2::new(
            (p_film_discrete.x + &self.filter_radius).ceil() + 1,
            (p_film_discrete.y + &self.filter_radius).ceil() + 1,
        );
        p0 = na::max(p0, &self.pixel_bounds.p_min as Point2<i32>);
        p1 = na::min(p1, &self.pixel_bounds.p_max as Point2<i32>);

        // precompute x and y filter table offsets
        // TODO is this correct?
        let ifx: Vec<i32> = (p0.x..p1.x)
            .map(|x| {
                let fx = ((x - p_film_discrete.x as i32)
                    * self.inv_filter_radius.x
                    * self.filter_table_size)
                    .abs();
                min(fx.floor() as i32, self.filter_table_size as i32)
            })
            .collect();
        let ify: Vec<i32> = (p0.y..p1.y)
            .map(|y| {
                let fy = ((y - p_film_discrete.y as i32)
                    * self.inv_filter_radius.y
                    * self.filter_table_size)
                    .abs();
                min(fy.floor() as i32, self.filter_table_size as i32)
            })
            .collect();

        // loop over filter support and add saple to pixel array
        for y in p0.y..p1.y {
            for x in p0.x..p1.x {
                // evaluate filter value at pixel (x, y)
                let offset = ify[y] * self.filter_table_size + ifx[x];
                let filter_weight = self.filter_table[offset];

                // update pixel values with filtered sample contribution
                let offset = self.get_pixel_offset(&Point2::new(x, y));
                let pixel = &mut self.pixels[offset as usize];
                pixel.contrib_sum += lighting * sample_weight * filter_weight;
                pixel.filter_weight_sum += filter_weight;
            }
        }
    }

    // Takes pixel coordinates w.r.t the image and converts them to coordinates in the film tile
    // before indexing into the pixels array
    pub fn get_pixel_offset(&self, p: &Point2<i32>) -> i32 {
        let width = self.pixel_bounds.p_max.x - self.pixel_bounds.p_min.x;
        let offset = (p.x - self.pixel_bounds.p_min.x) + (p.y - self.pixel_bounds.p_min.y) * width;

        offset
    }

    pub fn get_pixel_bounds(&self) -> Bounds2i {
        self.pixel_bounds
    }
}

// Models the sensing device in the simulated camera by finding the radiance for each camera ray
// and determines the sample's contribution to the pixels around the point on the film plane
// where the camera ray behan and updates its representation of the image
//
// Film is also responsible for determining the range of integer pixel values that the sampler
// is responsible for generating samples for..
pub struct Film {
    pub resolution: Point2<i32>,
    pub diagonal: Float,
    pub filter: Filter,
    pub filename: String,
    pub cropped_pixel_bounds: Bounds2i,
    pixels: Vec<Pixel>,
    filter_table: [Float; FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH], //wrong
    scale: Float,
}

impl Film {
    pub fn new(
        resolution: Point2<i32>,
        crop_window: Bounds2f,
        filter: Filter,
        diagonal: Float,
        filename: String,
        scale: Float,
    ) -> Film {
        // Compute film image bounds
        // TODO bounds should have been integer but couldn't get to work with nalgebra generics
        let cropped_pixel_bounds = Bounds2i::new(
            Point2::new(
                (resolution.x * crop_window.p_min.x).ceil(),
                (resolution.y * crop_window.p_min.y).ceil(),
            ),
            Point2::new(
                (resolution.x * crop_window.p_max.x).ceil(),
                (resolution.y * crop_window.p_max.y).ceil(),
            ),
        );

        //Allocate film image storage
        // TODO can do functionally instead of with mutation
        let pixels = Vec::with_capacity(cropped_pixel_bounds.area());

        // precompute filter weight table
        let mut filter_table = [0.0; FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH];
        let mut offset = 0;
        for y in 0..FILTER_TABLE_WIDTH {
            for x in 0..FILTER_TABLE_WIDTH {
                let p_x = (x as Float + 0.5) * filter.radius.x / FILTER_TABLE_WIDTH as Float;
                let p_y = (y as Float + 0.5) * filter.radius.y / FILTER_TABLE_WIDTH as Float;
                let p = Point2::new(p_x, p_y);
                filter_table[offset] = filter.evaluate(&p);
                offset += 1;
            }
        }

        Film {
            resolution,
            diagonal: diagonal * 0.001,
            filter,
            filename,
            cropped_pixel_bounds,
            pixels,
            filter_table,
            scale,
        }
    }

    // Returns the area to be sampled
    // Because pixel reconstruction spans a number of pixels, sampler will need to generate image
    // samples a bit outside of the range of pixels that will actually be otput
    pub fn get_sample_bounds(&self) -> Bounds2i {
        let p_min = Point2::new(
            (self.cropped_pixel_bounds.p_min.x + 0.5 - self.filter.radius).floor(),
            (self.cropped_pixel_bounds.p_min.y + 0.5 - self.filter.radius).floor(),
        );
        let p_max = Point2::new(
            (self.cropped_pixel_bounds.p_max.x - 0.5 + self.filter.radius).ceil(),
            (self.cropped_pixel_bounds.p_max.y - 0.5 + self.filter.radius).ceil(),
        );

        Bounds2i { p_min, p_max }
    }

    // Returns extent of the film in the scene
    pub fn get_physical_extent(&self) -> Bounds2f {
        let aspect = self.resolution.y / self.resolution.x;
        let x = ((self.diagonal * self.diagonal) / (1.0 + aspect * aspect)).sqrt();
        let y = aspect * x;

        Bounds2f::new(
            Point2::new(-x / 2.0, -y / 2.0),
            Point2::new(x / 2.0, y / 2.0),
        )
    }

    // TODO probably better to make Film Iterable and iterate over tile boundaries
    // Returns pointer to a FilmTile object that stores contributions for the pixels in the
    // corresponding region of the image. Ownership of the FilmTile and the data it stores
    // is exclusice to the caller
    pub fn get_film_tile(&self, sample_bounds: &Bounds2i) -> FilmTile {
        // bound image pixels that samples in sample_bounds contribute to
        let half_pixel = Vector2::new(0.5, 0.5);
        let float_bounds = Bounds2f::new(
            Point2::new(
                sample_bounds.p_min.x as Float,
                sample_bounds.p_min.y as Float,
            ),
            Point2::new(
                sample_bounds.p_max.x as Float,
                sample_bounds.p_max.y as Float,
            ),
        );
        let p0 = Point2::new(
            (float_bounds.p_min.x - half_pixel - self.filter.radius).ceil() as i32,
            (float_bounds.p_min.y - half_pixel - self.filter.radius).ceil() as i32,
        );
        let p1 = Point2::new(
            (float_bounds.p_max.x - half_pixel + self.filter.radius).floor() as i32 + 1,
            (float_bounds.p_max.y - half_pixel + self.filter.radius).floor() as i32 + 1,
        );
        let tile_pixel_bounds = self.cropped_pixel_bounds.intersect(&Bounds2i::new(p0, p1));

        FilmTile::New(
            tile_pixel_bounds,
            self.filter.radius,
            &self.filter_table,
            FILTER_TABLE_WIDTH,
        )
    }

    pub fn merge_film_tile(&mut self, tile: &FilmTile) {
        for pixel in tile.get_pixel_bounds() {
            // merge pixel into pixels
            let tile_offset = tile.get_pixel_offset(&pixel);
            let tile_pixel = tile.pixels[tile_offset as usize];
            let offset = self.get_pixel_offset(&pixel);
            let merge_pixel = &mut self.pixels[offset as usize];
            let xyz = tile_pixel.contrib_sum.to_xyz();
            for i in 0..3 {
                merge_pixel.xyz[i] += xyz[i];
            }
            merge_pixel.filter_weight_sum += tile_pixel.filter_weight_sum;
        }
    }

    pub fn set_image(&self, img: Vec<Spectrum>) {
        let n_pixels = self.cropped_pixel_bounds.surface_area();
        for i in 0..n_pixels {
            self.pixels[i] = Pixel::new(img[i], 1.0, [0.0; 3], 0.0);
        }
    }

    fn get_pixel_offset(&self, p: &Point2<i32>) -> i32 {
        let width = self.cropped_pixel_bounds.p_max.x - self.cropped_pixel_bounds.p_min.x;
        let offset = (p.x - self.cropped_pixel_bounds.p_min.x)
            + (p.y - self.cropped_pixel_bounds.p_min.y) * width;

        offset
    }
}
