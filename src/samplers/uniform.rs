use nalgebra::Point2;
use rand::prelude::*;

use crate::camera::CameraSample;
use crate::pbrt::Float;
use crate::random::Random;
use crate::sampler::{Sampler, SamplerIterator};

pub struct UniformSampler {
    pub samples_per_pixel: u64,
    rng: Random,
}

impl UniformSampler {
    pub fn new(samples_per_pixel: u64) -> UniformSampler {
        let rng = Random::new();
        UniformSampler {
            samples_per_pixel,
            rng,
        }
    }

    pub fn clone_with_seed(&self, seed: u64) -> UniformSampler {
        let mut sampler = UniformSampler::new(self.samples_per_pixel);
        sampler.rng.set_sequence(seed);
        sampler
    }

    pub fn reseed(&mut self, seed: u64) {
        self.rng.set_sequence(seed)
    }

    pub fn get_1d(&mut self) -> Float {
        self.rng.uniform_float()
    }

    pub fn get_2d(&mut self) -> Point2<Float> {
        Point2::new(self.rng.uniform_float(), self.rng.uniform_float())
    }

    pub fn sample(&mut self, pixel: &Point2<i32>) -> CameraSample {
        let sampled_2d = &self.get_2d();
        let p_film = Point2::new(
            pixel.x as Float + sampled_2d.x,
            pixel.y as Float + sampled_2d.y,
        );
        let time = self.get_1d();
        let p_lens = self.get_2d();
        let sample = CameraSample::new(p_film, p_lens, time);
        sample
    }
}
