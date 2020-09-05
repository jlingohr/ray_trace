extern crate nalgebra as na;

use na::geometry::{Point2, Point3};
use na::Vector3;

use crate::geometry::ray::Ray;
use crate::geometry::{face_forward, offset_ray_origin};
use crate::material::Material;
use crate::pbrt::{Float, Spectrum, SHADOW_EPSILON};

use self::na::{Matrix2, Vector2};

// fn face_forward<'a, 'b>(n: &'a Vector3<Float>, v: &'b Vector3<Float>) -> &'a Vector3<Float> {
//     if n.dot(&v) < 0.0 {
//         -n
//     } else {
//         n
//     }
// }

// Represent local information at a point on a 2D surface

// pub trait SurfaceInteraction {
//     fn isSurfaceInteraction(&self) -> bool;
//     fn spawnRay(&self, d: &Vector3<Float>) -> Ray;
//     fn spawnRayToPoint(&self, p2: &Point3<Float>) -> Ray;
//     fn spawnRayToInteraction(&self, it: &SurfaceInteraction) -> Ray;
//     fn isMediumInteraction(&self) -> bool;
//     fn getMedium(&self, w: &Vector3<Float>) -> Medium;
// }

#[derive(Debug, Copy, Clone)]
pub struct CommonInteraction {
    pub point: Point3<Float>,
    pub time: Float,
    pub p_error: Vector3<Float>,
    pub w0: Vector3<Float>,
    pub normal: Vector3<Float>,
}

impl CommonInteraction {
    pub fn new(
        point: Point3<Float>,
        time: Float,
        p_error: Vector3<Float>,
        w_out: Vector3<Float>,
        normal: Vector3<Float>,
    ) -> CommonInteraction {
        CommonInteraction {
            point,
            time,
            p_error,
            w0: w_out,
            normal,
        }
    }

    // pub fn set_face_forward(&mut self, n: &Vector3<Float>) {
    //     self.normal = *face_forward(self.normal, n);
    // }

    pub fn spawn_ray(&self, d: &Vector3<Float>) -> Ray {
        let o = offset_ray_origin(&self.point, &self.p_error, &self.normal, &d);
        Ray::new(o, d.clone_owned(), self.time, Float::INFINITY)
    }

    pub fn spawn_ray_to(&self, p2: &Point3<Float>) -> Ray {
        let w = p2 - &self.point;
        let origin = offset_ray_origin(&self.point, &self.p_error, &self.normal, &w);
        let d = p2 - &origin;
        Ray::new(origin, d, 1.0 - SHADOW_EPSILON, self.time)
    }

    pub fn spawn_ray_from_interaction(&self, it: &CommonInteraction) -> Ray {
        let w = &it.point - &self.point;
        let origin = offset_ray_origin(&self.point, &self.p_error, &self.normal, &w);
        let target = offset_ray_origin(&it.point, &it.p_error, &it.normal, &(&origin - &it.point));
        let d = target - &origin;
        Ray::new(origin, d, self.time, 1.0 - SHADOW_EPSILON)
    }
}

// Assumes for closed shapes that the normal is oriented such that it points to the
// outside of the shape
#[derive(Debug, Copy, Clone)]
pub struct Shading {
    pub normal: Vector3<Float>,
    pub dpdu: Vector3<Float>,
    pub dpdv: Vector3<Float>,
    pub dndu: Vector3<Float>,
    pub dndv: Vector3<Float>,
}

impl Shading {
    pub fn new(
        n: Vector3<Float>,
        dpdu: Vector3<Float>,
        dpdv: Vector3<Float>,
        dndu: Vector3<Float>,
        dndv: Vector3<Float>,
    ) -> Shading {
        Shading {
            normal: n,
            dpdu,
            dpdv,
            dndu,
            dndv,
        }
    }

    // fn set_face_forward(n: &Vector3<Float>) {
    //     self.n = face_forward(self.n, n);
    // }
}

// trait Shape {
//     fn reverse_orientation() -> bool;
//     fn transform_swap_handedness() -> bool;
// }

// Represents local information at a point on a 2D surface
pub struct SurfaceInteraction {
    pub interaction: CommonInteraction,
    // pub material: Material,
    pub uv: Point2<Float>,
    // (u, v) coords from parameterization of the surface
    pub dpdu: Vector3<Float>,
    pub dpdv: Vector3<Float>,
    // pub shape: &'a Shape,
    pub material: Material,
    pub dndu: Vector3<Float>,
    pub dndv: Vector3<Float>,
    pub shading: Shading,
    pub dpdx: Vector3<Float>,
    pub dpdy: Vector3<Float>,
    pub dudx: Float,
    pub dvdx: Float,
    pub dudy: Float,
    pub dvdy: Float,
    pub dist: Float,
}

impl SurfaceInteraction {
    pub fn new(
        p: Point3<Float>,
        p_error: Vector3<Float>,
        uv: Point2<Float>,
        wo: Vector3<Float>,
        dpdu: Vector3<Float>,
        dpdv: Vector3<Float>,
        dndu: Vector3<Float>,
        dndv: Vector3<Float>,
        time: Float,
        normal: Vector3<Float>,
        material: Material,
        dist: Float,
    ) -> SurfaceInteraction {
        let interaction = CommonInteraction::new(p, time, p_error, wo, normal.clone_owned()); //TODO is this ok?
                                                                                              //TODO better way to do this than cloning
        let shading = Shading::new(
            normal.clone_owned(),
            dpdu.clone_owned(),
            dpdv.clone_owned(),
            dndu.clone_owned(),
            dndv.clone_owned(),
        );
        SurfaceInteraction {
            interaction,
            // shape.get_material(),
            uv,
            dpdu,
            dpdv,
            // shape,
            material,
            dndu,
            dndv,
            shading,
            dpdx: Vector3::new(0.0, 0.0, 0.0),
            dpdy: Vector3::new(0.0, 0.0, 0.0),
            dudx: 0.0,
            dvdx: 0.0,
            dudy: 0.0,
            dvdy: 0.0,
            dist,
        }
    }

    // pub fn set_shading_geometry(
    //     &self,
    //     dpdus: &Vector3<Float>,
    //     dpdvs: &Vector3<Float>,
    //     dndus: &Vector3<Float>,
    //     dndvs: &Vector3<Float>,
    //     orientation_is_authoritative: bool,
    // ) {
    //     self.shading.normal = dpdus.cross(dpdvs).normalize();
    //     if self.shape.reverse_orientation ^ self.shape.transform_swap_handedness {
    //         self.shading.normal = -self.shading.normal;
    //     }
    //     if orientation_is_authoritative {
    //         self.interaction.set_face_forward(shading.n);
    //     } else {
    //         self.shading.set_face_wordward(self.interaction.n);
    //     }
    //     self.shading.dpdu = dpdus;
    //     self.shading.dpdv = dpdvs;
    //     self.shading.dndu = dndus;
    //     self.shading.dndv = dndvs;
    // }

    //
    pub fn compute_scattering_functions(&mut self, ray: &Ray) {
        self.compute_differentials(ray);
    }

    // compute emitted radiance at a surface point intersected by a ray
    pub fn le(&self, w: &Vector3<Float>) -> Spectrum {
        // let area: Option<AreaLight> = self.primitive.get_area_light();
        // area.map_or(Spectrum::new(0.0), |area| -> area.l(w))
        // todo need way to store primitive in interaction?
        Spectrum::new(0.0)
    }

    // Compute information about the projection size of the surface area around the intersection
    // on the image plane for use in texture antialiasing
    fn compute_differentials(&mut self, ray: &Ray) {
        match &ray.differential {
            Some(differential) => {
                // estimate screen space change in p and (u, v)
                let d = -self.interaction.normal.dot(&Vector3::new(
                    self.interaction.point.x,
                    self.interaction.point.y,
                    self.interaction.point.z,
                ));
                // TODO should be a proper way to make this converion
                let rx_origin_vec = Vector3::new(
                    differential.rx_origin.x,
                    differential.rx_origin.y,
                    differential.rx_origin.z,
                );
                let tx = (-self.interaction.normal.dot(&rx_origin_vec) / -d)
                    / self.interaction.normal.dot(&differential.rx_direction);
                let px = &differential.rx_origin + tx * &differential.rx_direction;

                let ry_origin_vec = Vector3::new(
                    differential.ry_origin.x,
                    differential.ry_origin.y,
                    differential.ry_origin.z,
                );
                let ty = (-self.interaction.normal.dot(&ry_origin_vec) - d)
                    / self.interaction.normal.dot(&differential.ry_direction);
                let py = &differential.ry_origin + tx * &differential.ry_direction;

                self.dpdx = px - &self.interaction.point;
                self.dpdy = py - &self.interaction.point;

                // compute (u, v) offsets at auxiliary point
                let dim = if self.interaction.normal.x.abs() > self.interaction.normal.y.abs()
                    && self.interaction.normal.x.abs() > self.interaction.normal.z.abs()
                {
                    [1, 2]
                } else if self.interaction.normal.y.abs() > self.interaction.normal.z.abs() {
                    [0, 2]
                } else {
                    [0, 1]
                };

                // initialize A, Bx, and By matrices for offset computation
                // TODO probably refactor this out into func
                let A = Matrix2::new(
                    self.dpdu[dim[0]],
                    self.dpdv[dim[0]],
                    self.dpdu[dim[1]],
                    self.dpdv[dim[1]],
                );
                let Bx = Vector2::new(
                    px[dim[0]] - &self.interaction.point[dim[0]],
                    px[dim[1]] - &self.interaction.point[dim[1]],
                );
                let By = Vector2::new(
                    py[dim[0]] - &self.interaction.point[dim[0]],
                    py[dim[1]] - &self.interaction.point[dim[1]],
                );

                let decomp = A.qr();
                let x: Option<Vector2<Float>> = decomp.solve(&Bx);
                let y: Option<Vector2<Float>> = decomp.solve(&By);

                match (x, y) {
                    (Some(x), Some(y)) => {
                        self.dudx = x.x;
                        self.dvdx = x.y;
                        self.dudy = y.x;
                        self.dvdy = y.y;
                    }

                    (Some(x), None) => {
                        self.dudx = x.x;
                        self.dvdx = x.y;
                        self.dudy = 0.0;
                        self.dvdy = 0.0;
                    }

                    (None, Some(y)) => {
                        self.dudx = 0.0;
                        self.dvdx = 0.0;
                        self.dudy = y.x;
                        self.dvdy = y.y;
                    }

                    (None, None) => {
                        self.dudx = 0.0;
                        self.dvdx = 0.0;
                        self.dudy = 0.0;
                        self.dvdy = 0.0;
                    }
                }
            }
            _ => {
                self.dpdx = Vector3::new(0.0, 0.0, 0.0);
                self.dpdy = Vector3::new(0.0, 0.0, 0.0);
                self.dudx = 0.0;
                self.dvdx = 0.0;
                self.dudy = 0.0;
                self.dvdy = 0.0;
            }
        }
    }
}
