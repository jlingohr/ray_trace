// // StratifiedSamplers will subdivide pixel areas into rectangular regions and generate a single
// // sample inside each region.
//
// use crate::pbrt::{Float, ONE_MINUS_EPSILON};
// use crate::random::Random;
// use nalgebra::Point2;
// use rand::seq::SliceRandom;
//
// fn stratified_sample_1d(i: usize, rng: &mut Random, jitter: bool) -> Float {
//     let inv_n_samples = 1.0 / n_samples;
//     let delta = if jitter { rng.uniform_float() } else { 0.5 };
//
//     std::cmp::min((i + delta) * inv_n_samples, ONE_MINUS_EPSILON)
// }
//
// fn stratified_sampled_2d(
//     x: usize,
//     y: usize,
//     nx: u32,
//     ny: u32,
//     rng: &mut Random,
//     jitter: bool,
// ) -> Point2<Float> {
//     let dx = 1.0 / nx as Float;
//     let dy = 1.0 / ny as Float;
//     let jx = if jitter { rng.uniform_float() } else { 0.5 };
//     let jy = if jitter { rng.uniform_float() } else { 0.5 };
//     let p = Point2::new(
//         std::cmp::min((x as Float + jx) * dx, ONE_MINUS_EPSILON),
//         std::cmp::min((y as Float + jy) * dy, ONE_MINUS_EPSILON),
//     );
//     p
// }
//
// struct StratifiedSampler {
//     samples_1d_array_sizes: Vec<i32>,
//     samples_2d_array_sizes: Vec<i32>,
//     sample_array_1d: Vec<Vec<Float>>,
//     sample_array_2d: Vec<Vec<Point2<Float>>>,
//     samples_x: u32,
//     samples_y: u32,
//     jitter_samples: bool,
// }
//
// impl StratifiedSampler {
//     pub fn new(
//         samples_x: u32,
//         samples_y: u32,
//         jitter_samples: bool,
//         n_sampled_dimensions: u32,
//     ) -> StratifiedSampler {
//         let samples_per_pixel = samples_x * samples_y;
//         let mut rng = Random::new();
//
//         let samples_1d_array_sizes: Vec<u32>; //TODO wtf is this
//         let sample_array_1d: Vec<Vec<Float>>;
//
//         // generate single stratified samples for the pixel
//         let samples_1d: Vec<Vec<Float>> = (0..n_sampled_dimensions)
//             .map(|x| {
//                 let mut v: Vec<Float> = (0..samples_per_pixel)
//                     .map(|y| stratified_sample_1d(y as usize, &mut rng.rng, jitter_samples))
//                     .collect();
//                 v.shuffle(&mut rng.rng)
//             })
//             .collect();
//
//         let samples_2d: Vec<Vec<Point2<Float>>> = (0..n_sampled_dimensions)
//             .map(|x| {
//                 let mut v: Vec<Point2<Float>> = (0..samples_per_pixel)
//                     .map(|y| {
//                         stratified_sampled_2d(
//                             x as usize,
//                             y as usize,
//                             samples_x,
//                             samples_y,
//                             &mut rng.rng,
//                             jitter_samples,
//                         )
//                     })
//                     .collect();
//                 v.shuffle(&mut rng.rng)
//             })
//             .collect();
//
//         //TODO this should be apart of StartPixel, or your iterator
//         // generate arrays of stratified samples for the pixel
//         for i in 0..samples_1d_array_sizes.len() {
//             for j in 0..samples_per_pixel {
//                 let count = samples_1d_array_sizes[i];
//                 sample_array_1d[i][(j * count) as usize]
//             }
//         }
//
//         StratifiedSampler {}
//     }
// }
