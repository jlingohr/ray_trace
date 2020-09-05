extern crate nalgebra as na;

use std::cmp::{max, min};

use na::geometry::Point3;
use na::Vector3;

use crate::efloat::{quadratic, EFloat};
use crate::geometry::bounds::Bounds3;
use crate::geometry::ray::Ray;
use crate::interaction::SurfaceInteraction;
use crate::material::Material;
use crate::pbrt::{clamp, distance, gamma, radians, Float, PI};
use crate::shape::Shape;
use crate::transforms::Transform;

use self::na::Point2;

// Reoresents a sphere centered at the origin in object space. To place elsewhere in
// the scene you must apply an appropriate transformation.
// TODO switch to RF instead of lifetime reference
#[derive(Debug, Copy, Clone)]
pub struct Sphere {
    pub object_to_world: Transform,
    pub world_to_object: Transform,
    pub reverse_orientation: bool,
    pub transform_swap_handedness: bool,
    pub material: Material,
    radius: Float,
    z_min: Float,
    z_max: Float,
    theta_min: Float,
    theta_max: Float,
    phi_max: Float,
}

impl Sphere {
    pub fn new(
        object_to_world: Transform,
        world_to_object: Transform,
        reverse_orientation: bool,
        radius: Float,
        z_min: Float,
        z_max: Float,
        phi_max: Float,
        material: Material,
    ) -> Sphere {
        // Truncate minimum and maximum z-values
        let z_min = clamp(z_min.min(z_max), -radius, radius);
        let z_max = clamp(z_min.max(z_max), -radius, radius);

        // Truncate sphereical coordinates
        let theta_min = clamp(z_min / radius, -1.0, 1.0).acos();
        let theta_max = clamp(z_max / radius, -1.0, 1.0).acos();
        let phi_max = radians(clamp(phi_max, 0.0, 360.0));

        Sphere {
            object_to_world,
            world_to_object,
            reverse_orientation,
            transform_swap_handedness: object_to_world.swap_handedness(),
            material,
            radius,
            z_min,
            z_max,
            theta_min,
            theta_max,
            phi_max,
        }
    }

    pub fn object_bound(self) -> Bounds3 {
        Bounds3::new(
            Point3::new(-self.radius, -self.radius, self.z_min),
            Point3::new(self.radius, self.radius, self.z_max),
        )
    }

    pub fn intersect(self, ray: &Ray) -> Option<SurfaceInteraction> {
        // transform ray to object space
        let (transformed_ray, o_error, d_error) =
            self.world_to_object.transform_ray_with_error(ray);

        // compute quadratic sphere coordinates
        // TODO use EFloat to accumulate errors
        let ox = EFloat::new(transformed_ray.origin.x as f32, o_error.x as f32);
        let oy = EFloat::new(transformed_ray.origin.y as f32, o_error.y as f32);
        let oz = EFloat::new(transformed_ray.origin.z as f32, o_error.z as f32);
        let dx = EFloat::new(transformed_ray.direction.x as f32, d_error.x as f32);
        let dy = EFloat::new(transformed_ray.direction.y as f32, d_error.y as f32);
        let dz = EFloat::new(transformed_ray.direction.z as f32, d_error.z as f32);

        let a = dx * dx + dy * dy + dz * dz;
        let b: EFloat = (dx * ox + dy * oy + dz * oz) * 2.0;
        let c = ox * ox + oy * oy + oz * oz
            - EFloat::new(self.radius as f32, 0.0) * EFloat::new(self.radius as f32, 0.0);

        // solve quadrative equation for t values
        match quadratic(a, b, c) {
            //TODO a lot of the code could be refactored since a lot of repetition
            Some((t0, t1)) => {
                // check quadratic shape t0 and t1 for nearest intersection
                if t0.upper_bound() > transformed_ray.t_max as f32 || t1.lower_bound() <= 0.0 {
                    return None;
                } else {
                    let mut t_shape_hit = t0;
                    if t_shape_hit.lower_bound() <= 0.0 {
                        t_shape_hit = t1;
                        if t_shape_hit.upper_bound() > ray.t_max as f32 {
                            return None;
                        }
                    }
                    let mut p_hit = transformed_ray.at(t_shape_hit.v as Float);

                    // refine sphere intersection point
                    p_hit = &p_hit * (self.radius / distance(&p_hit, &Point3::new(0.0, 0.0, 0.0)));
                    if p_hit.x == 0.0 && p_hit.y == 0.0 {
                        p_hit.x = 1e-5 * self.radius;
                    }
                    let mut phi = p_hit.y.atan2(p_hit.x);
                    if phi < 0.0 {
                        phi += 2.0 * PI;
                    }

                    // test sphere intersection against clipping parameters
                    if (self.z_min > -self.radius && p_hit.z < self.z_min)
                        || (self.z_max < self.radius && p_hit.z > self.z_max)
                        || phi > self.phi_max
                    {
                        if t_shape_hit == t1 {
                            return None;
                        }
                        if t1.upper_bound() > transformed_ray.t_max as f32 {
                            return None;
                        }
                        t_shape_hit = t1;
                        p_hit = transformed_ray.at(t_shape_hit.v as Float);
                        p_hit *= self.radius / distance(&p_hit, &Point3::new(0.0, 0.0, 0.0));
                        if p_hit.x == 0.0 && p_hit.y == 0.0 {
                            p_hit.x = 1e-5 * self.radius;
                        }
                        phi = p_hit.y.atan2(p_hit.x);
                        if phi < 0.0 {
                            phi += 2.0 * PI;
                        }
                        if (self.z_min > -self.radius && p_hit.z < self.z_min)
                            || (self.z_min < self.radius && p_hit.z > self.z_max)
                            || phi > self.z_max
                        {
                            return None;
                        }
                    }
                    // find parametric representation of sphere hit
                    let u = phi / self.phi_max;
                    let theta = clamp(p_hit.z / self.radius, -1.0, 1.0).acos();
                    let v = (theta - self.theta_min) / (self.theta_max - self.theta_min);

                    // compute sphere dp/du and dp/dv
                    let z_radius = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
                    let inv_z_radius = 1.0 / z_radius;
                    let cos_phi = p_hit.x * inv_z_radius;
                    let sin_phi = p_hit.y * inv_z_radius;
                    let dpdu = Vector3::new(-self.phi_max * p_hit.y, self.phi_max * p_hit.x, 0.0);
                    let dpdv = (self.theta_max - self.theta_min)
                        * Vector3::new(
                            p_hit.z * cos_phi,
                            p_hit.z * sin_phi,
                            -self.radius * theta.sin(),
                        );

                    // compute sphere dn/du and dn/dv
                    let d2Pduu = -self.phi_max * self.phi_max * Vector3::new(p_hit.x, p_hit.y, 0.0);
                    let d2Pduv = (self.theta_max - self.theta_min)
                        * p_hit.z
                        * self.phi_max
                        * Vector3::new(-sin_phi, cos_phi, 0.0);
                    let d2Pdvv = -(self.theta_max - self.theta_min)
                        * (self.theta_max - self.theta_min)
                        * Vector3::new(p_hit.x, p_hit.y, p_hit.z);

                    // compute coefficients for fundamental forms
                    let E = dpdu.dot(&dpdu);
                    let F = dpdu.dot(&dpdv);
                    let G = dpdv.dot(&dpdv);
                    let normal = dpdu.cross(&dpdv).normalize();
                    let e = normal.dot(&d2Pduu);
                    let f = normal.dot(&d2Pduv);
                    let g = normal.dot(&d2Pdvv);

                    let inv_EGF2 = 1.0 / (E * G - F * F);
                    let dndu =
                        (f * f - e * G) * inv_EGF2 * dpdu + (e * F - f * E) * inv_EGF2 * dpdv;
                    let dndv =
                        (g * F - f * G) * inv_EGF2 * dpdu + (f * F - g * E) * inv_EGF2 * dpdv;

                    // copute error bounds for sphere intersection
                    let p_error = gamma(5) as Float
                        * Vector3::new(p_hit.x.abs(), p_hit.y.abs(), p_hit.z.abs());

                    // initialize surface interaction from parametric information
                    // update hit for quadratic intersection
                    //TODO also want to return t_hit to allow subsequent intersection tests to terminate early
                    // if the potential hit would be farther away than the existing intersection
                    let t_hit = t_shape_hit.v as Float;
                    let mut n = dpdu.cross(&dpdv).normalize();
                    // adjust normal based on orientation and handedness
                    if self.reverse_orientation ^ self.transform_swap_handedness {
                        n *= -1.0;
                    }
                    let mut isect = SurfaceInteraction::new(
                        p_hit,
                        p_error,
                        Point2::new(u, v),
                        -transformed_ray.direction,
                        dpdu,
                        dpdv,
                        dndu,
                        dndv,
                        transformed_ray.time,
                        n,
                        self.material,
                        t_hit,
                    );

                    // convert to world space
                    self.object_to_world
                        .transform_surface_interaction(&mut isect);
                    Some(isect)
                }
            }
            None => None,
        }
    }

    // similar to intersect but does not initialize a SurfaceInteraction
    // TODO do we really need to recompute all this crap?
    pub fn intersect_p(self, ray: &Ray) -> bool {
        // transform ray to object space
        let (transformed_ray, o_error, d_error) =
            self.world_to_object.transform_ray_with_error(ray);

        // compute quadratic sphere coordinates
        // TODO use EFloat to accumulate errors
        let ox = EFloat::new(transformed_ray.origin.x as f32, o_error.x as f32);
        let oy = EFloat::new(transformed_ray.origin.y as f32, o_error.y as f32);
        let oz = EFloat::new(transformed_ray.origin.z as f32, o_error.z as f32);
        let dx = EFloat::new(transformed_ray.direction.x as f32, d_error.x as f32);
        let dy = EFloat::new(transformed_ray.direction.y as f32, d_error.y as f32);
        let dz = EFloat::new(transformed_ray.direction.z as f32, d_error.z as f32);

        let a = dx * dx + dy * dy + dz * dz;
        let b = (dx * ox + dy * oy + dz * oz) * 2.0;
        let c = ox * ox + oy * oy + oz * oz
            - EFloat::new(self.radius as f32, 0.0) * EFloat::new(self.radius as f32, 0.0);

        // solve quadrative equation for t values
        match quadratic(a, b, c) {
            //TODO a lot of the code could be refactored since a lot of repetition
            Some((t0, t1)) => {
                // check quadratic shape t0 and t1 for nearest intersection
                if t0.upper_bound() > transformed_ray.t_max as f32 || t1.lower_bound() <= 0.0 {
                    return false;
                } else {
                    let mut t_shape_hit = t0;
                    if t_shape_hit.lower_bound() <= 0.0 {
                        t_shape_hit = t1;
                        if t_shape_hit.upper_bound() > ray.t_max as f32 {
                            return false;
                        }
                    }
                    let mut p_hit = transformed_ray.at(t_shape_hit.v as Float);

                    // refine sphere intersection point
                    p_hit = &p_hit * (self.radius / distance(&p_hit, &Point3::new(0.0, 0.0, 0.0)));
                    if p_hit.x == 0.0 && p_hit.y == 0.0 {
                        p_hit.x = 1e-5 * self.radius;
                    }
                    let mut phi = p_hit.y.atan2(p_hit.x);
                    if phi < 0.0 {
                        phi += 2.0 * PI;
                    }

                    // test sphere intersection against clipping parameters
                    if (self.z_min > -self.radius && p_hit.z < self.z_min)
                        || (self.z_max < self.radius && p_hit.z > self.z_max)
                        || phi > self.phi_max
                    {
                        if t_shape_hit == t1 {
                            return false;
                        }
                        if t1.upper_bound() > transformed_ray.t_max as f32 {
                            return false;
                        }
                        t_shape_hit = t1;
                        p_hit = transformed_ray.at(t_shape_hit.v as Float);
                        p_hit *= self.radius / distance(&p_hit, &Point3::new(0.0, 0.0, 0.0));
                        if p_hit.x == 0.0 && p_hit.y == 0.0 {
                            p_hit.x = 1e-5 * self.radius;
                        }
                        phi = p_hit.y.atan2(p_hit.x);
                        if phi < 0.0 {
                            phi += 2.0 * PI;
                        }
                        if (self.z_min > -self.radius && p_hit.z < self.z_min)
                            || (self.z_min < self.radius && p_hit.z > self.z_max)
                            || phi > self.z_max
                        {
                            return false;
                        }
                    }

                    true
                }
            }
            None => false,
        }
    }

    pub fn area(&self) -> Float {
        self.phi_max * self.radius * (self.z_max - self.z_min)
    }
}
