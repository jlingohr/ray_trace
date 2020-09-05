use nalgebra::Point2;

use crate::camera::CameraSample;
use crate::film::Pixel;
use crate::pbrt::Float;
use crate::random::Random;
use crate::samplers::uniform::UniformSampler;

enum SamplerType {}

struct SamplerParams {
    num_samples: i64,
    sampler_type: SamplerType,
}

// impl SamplerParams {
//     pub fn to_sampler(&self, p: &Point2<i32>) -> Sampler {
//         unimplemented!()
//     }
// }

// Sampler is responsible for
// 1) choosing the points on the image plane from which rays are traced
// 2) supplying the sample positions used by the integrators for estimating the value of
//      the light transport integral
//
// If arrays of samples are needed, they must be requested before rendering begins

pub trait Sampler {
    // TODO probably better to just initialize a new sampler with these parameters
    // fn start_pixel(&self, p: &Point2<u32>) {
    //     self.current_pixel = p;
    //     self.current_pixel_sample_index = 0;
    //     self.array_1d_offset = 0;
    //     self.array_2d_offset = 0;
    // }

    // When work for one sample is complete, the integrator calls start_next_sample()
    // which notifes the Sampler that subsequent request for sample components should
    // return values starting at the first dimension of the next sample for the current pixel. Should
    // return true until number of originally request samples per pixel has been generated
    // TODO probably just use Iterator?

    // fn start_next_sample(&self) -> bool {
    //     self.array_1d_offset = 0;
    //     self.array_2y_offset = 0;
    //     self.current_pixel_sample_index += 1;
    //
    //     self.current_pixel_sample_index < self.samples_per_pixel
    // }
    //
    // fn set_sample_num(&self, sample_num: u64) -> bool {
    //     self.array_1d_offset = 0;
    //     self.array_2y_offset = 0;
    //     self.current_pixel_sample_index = sample_num;
    //
    //     self.current_pixel_sample_index < self.samples_per_pixel
    // }

    // Returns the sample value for the next dimensio of the current sapmle vector
    fn get_1d(&self) -> Float;

    // Returns the sample values
    fn get_2d(&self) -> Point2<Float>;

    // initializes a CameraSample for a given pixel
    // fn get_camera_sample(&self, p_raster: &Point2<u32>) -> CameraSample {
    //     let p_film = p_raster + self.get_2d();
    //     let time = self.get_1d();
    //     let p_lens = self.get_2d();
    //     CameraSample::new(p_film, p_lens, time)
    // }

    fn request_1d_array(&self, n: u32);

    fn request_2d_array(&self, n: u32);

    // Returns array of samples whose size is given by n
    fn get_1d_array(&self, n: u32) -> Vec<Float>;

    fn get_2d_array(&self, n: u32) -> Vec<Point2<Float>>;

    // Reset current_[12]d_dim
    fn reset(&self);
}

pub struct SamplerIterator<'a> {
    pixel: Point2<Float>,
    samples_per_pixel: u64,
    // current_pixel: Point2<u32>,
    current_pixel_sample_index: u64,
    // array_1d_offset: usize,
    // array_2d_offset: usize,
    sampler: &'a mut UniformSampler,
}

impl<'a> SamplerIterator<'a> {
    pub fn new(
        pixel: Point2<Float>,
        samples_per_pixel: u64,
        sampler: &'a mut UniformSampler,
    ) -> SamplerIterator<'a> {
        SamplerIterator {
            pixel,
            samples_per_pixel,
            current_pixel_sample_index: 0,
            sampler,
        }
    }
}

impl<'a> Iterator for SamplerIterator<'a> {
    type Item = CameraSample;

    fn next(&mut self) -> Option<Self::Item> {
        // self.sampler.reset();
        // self.array_1d_offset = 0;
        // self.array_2d_offset = 0;
        // self.current_pixel_sample_index += 1;

        if self.current_pixel_sample_index <= self.samples_per_pixel {
            let sampled_2d = self.sampler.get_2d();
            let p_film = Point2::new(&self.pixel.x + sampled_2d.x, &self.pixel.y + sampled_2d.y);
            let time = self.sampler.get_1d();
            let p_lens = self.sampler.get_2d();
            let sample = CameraSample::new(p_film, p_lens, time);
            Some(sample)
        } else {
            None
        }
    }
}

// pub struct PixelSampler {
//     samples_per_pixel: u64,
//     samples_1d: Vec<Vec<Float>>,
//     samples_2d: Vec<Vec<Point2<Float>>>,
//     current_1d_dim: usize,
//     current_2d_dim: usize,
//     rng: Random,
// }
//
// impl PixelSampler {
//     pub fn new(samples_per_pixel: u64, n_sampled_dimensions: i32) -> PixelSampler {
//         let mut samples_1d: Vec<Vec<Float>> = Vec::new();
//         let mut samples_2d: Vec<Vec<Point2<Float>>> = Vec::new();
//         for i in 0..n_sampled_dimensions {
//             samples_1d.push(Vec::with_capacity(samples_per_pixel as usize));
//             samples_2d.push(Vec::with_capacity(samples_per_pixel as usize));
//         }
//         let rng = Random::new();
//
//
//
//         PixelSampler {
//             samples_per_pixel,
//             samples_1d,
//             samples_2d,
//             current_1d_dim: 0,
//             current_2d_dim: 0,
//             rng,
//         }
//     }
//
//     fn get_1d(&mut self, current_pixel_sample_index: usize) -> Float {
//         if self.current_1d_dim < self.samples_1d.len() {
//             let res = self.samples_1d[self.current_1d_dim][current_pixel_sample_index];
//             self.current_1d_dim += 1;
//             res
//         } else {
//             self.rng.uniform_float()
//         }
//     }
//
//     fn get_2d(&mut self, current_pixel_sample_index: Float) -> Point2<Float> {
//         if self.current_2d_dim < self.samples_2d.len() {
//             let res = self.samples_2d[self.current_2d_dim][current_pixel_sample_index];
//             self.current_2d_dim += 1;
//             res
//         } else {
//             Point2::new(self.rng.uniform_float(), self.rng_uniform_float()) //TODO
//         }
//     }
//
//     fn reset(&mut self) {
//         self.current_1d_dim = 0;
//         self.current_2d_dim = 0;
//     }
//
//
// }
