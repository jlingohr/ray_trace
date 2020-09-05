extern crate nalgebra as na;

use std::sync::Arc;

use na::geometry::{Point2, Point3};
use na::Vector3;

use crate::film::Film;
use crate::geometry::bounds::Bounds2f;
use crate::geometry::ray::Ray;
use crate::pbrt::{lerp, Float};
use crate::sampling::concentric_sample_disk;
use crate::transforms::{AnimatedTransform, Transform};

// use crate::interaction::Interaction;

// Hold all sample values needed to specify a camera ray
// Fields:
// * pfilm: point on film on which the generated ray carries radiance
// * plens: point on lens on through which the ray passes
// * time: time at which the ray should sample the scene
pub struct CameraSample {
    pub p_film: Point2<Float>,
    pub p_lens: Point2<Float>,
    pub time: Float,
}

impl CameraSample {
    pub fn new(pfilm: Point2<Float>, plens: Point2<Float>, time: Float) -> CameraSample {
        CameraSample {
            p_film: pfilm,
            p_lens: plens,
            time,
        }
    }
}

pub struct BaseCamera {
    camera_to_world: AnimatedTransform,
    shutter_open: Float,
    shutter_close: Float,
    // film: Arc<Film>,
    // medium: Medium,
}

impl BaseCamera {
    pub fn new(
        camera_to_world: AnimatedTransform,
        shutter_open: Float,
        shutter_close: Float, // film: Arc<Film>,
    ) -> BaseCamera {
        BaseCamera {
            camera_to_world,
            shutter_open,
            shutter_close,
            // film,
        }
    }
}

pub enum Camera {
    OrthographicCamera(OrthographicCamera),
}

impl Camera {
    pub fn generate_ray(&self, sample: CameraSample) -> Ray {
        match self {
            Camera::OrthographicCamera(camera) => camera.generate_ray(&sample),
        }
    }

    // Generates main camera way
    pub fn generate_ray_differential(&self, sample: &CameraSample) -> (Ray, Float) {
        match self {
            Camera::OrthographicCamera(camera) => camera.generate_ray_differential(sample),
        }
    }

    pub fn we(&self, ray: &Ray, p_raster2: Point2<Float>) {
        unimplemented!()
    }
    pub fn pdf_we(&self, ray: &Ray, pdf_pos: Float, pdf_dir: Float) {
        unimplemented!()
    }
    // fn sample_wi<T: Interaction>(&self, reference: &T, u: Point2<Float>, wi: Vector3<Float>, pdf: Float, p_raster: Point2<Float>, vis: VisibilityTester);

    pub fn get_film(&self) -> Arc<Film> {
        match self {
            Camera::OrthographicCamera(camera) => camera.get_film(),
        }
    }

    pub fn save_image(&self) {
        match self {
            Camera::OrthographicCamera(camera) => camera.save_image(),
        }
    }
}

pub struct OrthographicCamera {
    pub base_camera: BaseCamera,
    pub camera_to_screen: Transform,
    pub raster_to_camera: Transform,
    pub screen_to_raster: Transform,
    pub raster_to_screen: Transform,
    pub lens_radius: Float,
    pub focal_distance: Float,
    film: Arc<Film>,
    dx: Vector3<Float>,
    dy: Vector3<Float>,
}

impl OrthographicCamera {
    pub fn new(
        camera_to_world: AnimatedTransform,
        camera_to_screen: Transform,
        screen_window: Bounds2f,
        shutter_open: Float,
        shutter_close: Float,
        lens_radius: Float,
        focal_distance: Float,
        film: Arc<Film>,
    ) -> OrthographicCamera {
        let camera = BaseCamera::new(camera_to_world, shutter_open, shutter_close);
        let screen_to_raster =
            Transform::scale(film.resolution.x as Float, film.resolution.y as Float, 1.0)
                * Transform::scale(
                    1.0 / (screen_window.p_max.x - screen_window.p_min.x),
                    1.0 / (screen_window.p_min.y - screen_window.p_max.y),
                    1.0,
                )
                * Transform::translate(&Vector3::new(
                    -screen_window.p_min.x,
                    -screen_window.p_max.y,
                    0.0,
                ));
        let raster_to_screen = screen_to_raster.inverse();
        let raster_to_camera = camera_to_screen.inverse() * raster_to_screen;
        let dx = raster_to_camera.transform_vector(&Vector3::new(1.0, 0.0, 0.0));
        let dy = raster_to_camera.transform_vector(&Vector3::new(0.0, 1.0, 0.0));

        OrthographicCamera {
            base_camera: camera,
            camera_to_screen,
            raster_to_camera,
            screen_to_raster,
            raster_to_screen,
            lens_radius,
            focal_distance,
            film,
            dx,
            dy,
        }
    }

    fn generate_ray(&self, sample: &CameraSample) -> Ray {
        let p_film = Point3::new(sample.p_film.x, sample.p_film.y, 0.0);
        let p_camera = self.raster_to_camera.transform_point(&p_film);
        let time = lerp(
            sample.time,
            self.base_camera.shutter_open,
            self.base_camera.shutter_close,
        );

        // Modify ray for depth of field
        let ray = Ray::new(p_camera, Vector3::new(0.0, 0.0, 1.0), time, Float::INFINITY);
        if self.lens_radius > 0.0 {
            // sample point on lens
            let p_lens = self.lens_radius * concentric_sample_disk(&sample.p_lens);
            // compute point on plane of focus
            let ft = self.focal_distance / ray.direction.z;
            let p_focus = ray.at(ft);
            // update ray for effect of lens
            let origin = Point3::new(p_lens.x, p_lens.y, 0.0);
            let direction = (p_focus - ray.origin).normalize();

            let new_ray = Ray::new(origin, direction, time, ray.t_max);
            self.base_camera.camera_to_world.transform_ray(&new_ray)
        } else {
            let origin = self.raster_to_camera.transform_point(&p_film);
            let new_ray = Ray::new(origin, Vector3::new(0.0, 0.0, 1.0), time, ray.t_max);
            self.base_camera.camera_to_world.transform_ray(&new_ray)
        }
    }

    // TODO probably return the Ray instead
    fn generate_ray_differential(&self, sample: &CameraSample) -> (Ray, Float) {
        // compute raster and camera sample positions
        let p_film = Point3::new(sample.p_film.x, sample.p_film.y, 0.0);
        let p_camera = self.raster_to_camera.transform_point(&p_film);
        let time = lerp(
            sample.time,
            self.base_camera.shutter_open,
            self.base_camera.shutter_close,
        );

        // Modify ray for depth of field
        let mut ray = Ray::new(p_camera, Vector3::new(0.0, 0.0, 1.0), time, Float::INFINITY);
        ray = if self.lens_radius > 0.0 {
            // sample point on lens
            let p_lens = self.lens_radius * concentric_sample_disk(&sample.p_lens);
            // compute point on plane of focus
            let ft = self.focal_distance / ray.direction.z;
            let p_focus = ray.at(ft);
            // update ray for effect of lens
            let origin = Point3::new(p_lens.x, p_lens.y, 0.0);
            let direction = (p_focus - ray.origin).normalize();

            Ray::new(origin, direction, time, ray.t_max)
        } else {
            let origin = self.raster_to_camera.transform_point(&p_film);
            Ray::new(origin, Vector3::new(0.0, 0.0, 1.0), time, ray.t_max)
        };

        // compute ray differentials
        // compute differential changes in origin
        let dx_camera = self
            .raster_to_camera
            .transform_vector(&Vector3::new(1.0, 0.0, 0.0));
        let dy_camera = self
            .raster_to_camera
            .transform_vector(&Vector3::new(0.0, 1.0, 0.0));

        if self.lens_radius > 0.0 {
            let p_lens = self.lens_radius * concentric_sample_disk(&sample.p_lens);
            // compute point on plane of focus
            let ft = self.focal_distance / ray.direction.z;
            let p_focus = &p_camera + dx_camera + (ft * Vector3::new(0.0, 0.0, 1.0));
            let rx_origin = Point3::new(p_lens.x, p_lens.y, 0.0);
            let rx_direction = (p_focus - &rx_origin).normalize();

            let p_focus = &p_camera + dy_camera + (ft * Vector3::new(0.0, 0.0, 1.0));
            let ry_origin = Point3::new(p_lens.x, p_lens.y, 0.0);
            let ry_direction = (p_focus - &ry_origin).normalize();

            let new_ray = ray.with_differential(rx_origin, ry_origin, rx_direction, ry_direction);
            (
                self.base_camera.camera_to_world.transform_ray(&new_ray),
                1.0,
            )
        } else {
            let rx_origin = &ray.origin + dx_camera;
            let ry_origin = &ray.origin + dy_camera;
            let rx_direction = ray.direction.clone_owned();
            let ry_direction = ray.direction.clone_owned();

            let new_ray = ray.with_differential(rx_origin, ry_origin, rx_direction, ry_direction);
            (
                self.base_camera.camera_to_world.transform_ray(&new_ray),
                1.0,
            )
        }
    }

    fn we(&self, ray: &Ray, p_raster2: Point2<Float>) {
        unimplemented!()
    }

    fn pdf_we(&self, ray: &Ray, pdf_pos: Float, pdf_dir: Float) {
        unimplemented!()
    }

    // fn sample_wi<T: Interaction>(&self, reference: &T, u: Point2<Float>, wi: Vector3<Float>, pdf: Float, p_raster: Point2<Float>, vis: _) {
    //     unimplemented!()
    // }

    fn get_film(&self) -> Arc<Film> {
        self.film.clone()
    }

    pub fn save_image(&self) {
        //TODO where to get splat_scale?
        self.film.save_image(1.0);
    }

    // pub fn merge_film_tile(&mut self, film_tile: &FilmTile) {
    //     let mut film = self.base_camera.film;
    //     film.merge_film_tile(film_tile);
    // }
}
