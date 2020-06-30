use super::point::Point;
use super::vector::Vec3;

use rand::Rng;

const POINT_COUNT: usize = 256;

fn perlin_generate_perm() -> Vec<usize> {
    let mut p: Vec<usize> = Vec::with_capacity(POINT_COUNT);
    for i in 0..POINT_COUNT {
        p.push(i)
    }

    permute(&mut p);
    p
}

fn permute(p: &mut Vec<usize>) {
    let mut rng = rand::thread_rng();
    for i in (1..p.len()).rev() {
        p.swap(i, rng.gen_range(0, i));
    }
}

fn trilinear_interp(c: &[[[Vec3; 2]; 2]; 2], u: f64, v: f64, w: f64) -> f64 {
    let uu = u * u * (3.0 - (2.0 * u));
    let vv = v * v * (3.0 - (2.0 * v));
    let ww = w * w * (3.0 - (2.0 * w));
    let mut acc = 0.0;
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                let weight = Vec3::new(u - i as f64, v - j as f64, w - k as f64);
                acc += ((i as f64 * uu) + ((1 - i) as f64 * (1.0 - uu)))
                    * ((j as f64 * vv) + ((1 - j) as f64 * (1.0 - vv)))
                    * ((k as f64 * ww) + ((1 - k) as f64 * (1.0 - ww)))
                    * c[i][j][k].dot(&weight);
            }
        }
    }
    acc
}

#[derive(Clone)]
pub struct PerlinNoise {
    rand_vec: Vec<Vec3>,
    perm_x: Vec<usize>,
    perm_y: Vec<usize>,
    perm_z: Vec<usize>,
}

impl PerlinNoise {
    pub fn new() -> PerlinNoise {
        let mut rand_vec = Vec::with_capacity(POINT_COUNT);
        let mut rng = rand::thread_rng();
        for _ in 0..POINT_COUNT {
            rand_vec.push(Vec3::random_in_range(&mut rng, -1.0, 1.0).unit());
        }

        let perm_x = perlin_generate_perm();
        let perm_y = perlin_generate_perm();
        let perm_z = perlin_generate_perm();

        PerlinNoise {
            rand_vec,
            perm_x,
            perm_y,
            perm_z,
        }
    }

    pub fn noise(&self, p: &Point) -> f64 {
        let u = p.x - p.x.floor();
        let v = p.y - p.y.floor();
        let w = p.z - p.z.floor();

        let i = p.x.floor();
        let j = p.y.floor();
        let k = p.z.floor();

        let mut c: [[[Vec3; 2]; 2]; 2] = [[[Vec3::new(0.0, 0.0, 0.0); 2]; 2]; 2];

        for di in 0..2 {
            for dj in 0..2 {
                for dk in 0..2 {
                    let ix = self.perm_x[((i as i64 + di as i64) & 255) as usize];
                    let iy = self.perm_y[((j as i64 + dj as i64) & 255) as usize];
                    let iz = self.perm_z[((k as i64 + dk as i64) & 255) as usize];
                    c[di][dj][dk] = self.rand_vec[ix ^ iy ^ iz];
                }
            }
        }

        trilinear_interp(&c, u, v, w)
    }

    pub fn turb(&self, p: &Point, depth: usize) -> f64 {
        let mut acc = 0.0;
        let mut temp_p = *p;
        let mut weight = 1.0;

        for i in 0..depth {
            acc += weight * self.noise(&temp_p);
            weight *= 0.5;
            temp_p = 2.0 * (*p);
        }
        f64::abs(acc)
    }
}
