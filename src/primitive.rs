use crate::geometry::bounds::Bounds3;
use crate::geometry::ray::Ray;
use crate::interaction::SurfaceInteraction;
use crate::light::AreaLight;
use crate::material::Material;
use crate::pbrt::Float;
use crate::reflection::BSDF;
use crate::shape::Shape;

// pub trait Primitive: Intersectable + HasLight + HasMaterial {}

pub enum Primitive {
    Geometric(GeometricPrimitive),
}

impl Primitive {
    // Returns referencee to object describing the primitive's emission distribution, if the
    // the primitive is itself a light source. If the primitive is not emissive returns None
    fn get_area_light(&self) -> Option<AreaLight> {
        match self {
            Primitive::Geometric(primitive) => primitive.get_area_light(),
        }
    }

    // Returns material assigned to the primitive. Serves only to delineate a volume of
    // space for participating media. Should also be used to check if two rays have
    // intersected the same object by comparing their material pointers.
    fn get_material(&self) -> Material {
        match self {
            Primitive::Geometric(primitive) => primitive.get_material(),
        }
    }

    // Initializes representations of the light-scattering properties of the material at the
    // intersection point on the surface
    fn compute_scattering_functions(
        &self,
        isect: &SurfaceInteraction,
        allow_multiple_lobes: bool,
    ) -> Option<BSDF> {
        match self {
            Primitive::Geometric(primitive) => {
                primitive.compute_scattering_functions(isect, allow_multiple_lobes)
            }
        }
    }

    // Returns a box that encloses the primitive's geometry i world space
    fn world_bounds(&self) -> Bounds3 {
        match self {
            Primitive::Geometric(primitive) => self.world_bounds(),
        }
    }

    // ray intersection test
    fn intersect(&self, r: &Ray) -> Option<SurfaceInteraction> {
        match self {
            Primitive::Geometric(primitive) => primitive.intersect(r),
        }
    }

    fn intersect_p(&self, r: &Ray) -> bool {
        match self {
            Primitive::Geometric(primitive) => primitive.intersect_p(r),
        }
    }
}

pub struct GeometricPrimitive {
    shape: Shape,
    material: Material,
    area_light: Option<AreaLight>,
    // medium: Medium,
}

impl GeometricPrimitive {
    pub fn new(
        shape: Shape,
        material: Material,
        area_light: Option<AreaLight>,
    ) -> GeometricPrimitive {
        GeometricPrimitive {
            shape,
            material,
            area_light,
        }
    }

    pub fn world_bounds(&self) -> Bounds3 {
        unimplemented!()
    }

    pub fn intersect(&self, r: &Ray) -> Option<SurfaceInteraction> {
        self.shape.intersect(r)
        // if let Some(isect) = self.shape.intersect(r) {
        //     r.t_max = isect.interaction.dist;
        //     isect.primitive = self;
        //     if self.medium_interface.is_medium_transition() {
        //         isect.medium_interface = self.medium_interface;
        //     } else {
        //         isect.medium_interface = MediumInterface::new(r.medium);
        //     }
        // } else {
        //     None
        // }
    }

    pub fn intersect_p(&self, r: &Ray) -> bool {
        self.shape.intersect_p(r)
    }

    pub fn get_area_light(&self) -> Option<AreaLight> {
        unimplemented!()
    }

    pub fn get_material(&self) -> Material {
        unimplemented!()
    }

    pub fn compute_scattering_functions(
        &self,
        isect: &SurfaceInteraction,
        allow_multiple_lobes: bool,
    ) -> Option<BSDF> {
        self.material
            .compute_scattering_functions(isect, allow_multiple_lobes)
    }
}

// Groups multiple primitive objects together
// these should implement Primitive but never call get_area_light, get_material,
// or compute_scattering_functions
pub enum Aggregate {
    HittableList(HittableList),
}

impl Aggregate {
    pub fn world_bound(&self) -> Bounds3 {
        match self {
            Aggregate::HittableList(agg) => agg.world_bounds(),
        }
    }

    pub fn intersect(&self, r: &Ray) -> Option<SurfaceInteraction> {
        match self {
            Aggregate::HittableList(agg) => agg.intersect(r),
        }
    }

    pub fn intersect_p(&self, r: &Ray) -> bool {
        match self {
            Aggregate::HittableList(agg) => agg.intersect_p(r),
        }
    }
}

pub struct HittableList {
    objects: Vec<Primitive>,
}

impl HittableList {
    pub fn new(objects: Vec<Primitive>) -> HittableList {
        HittableList { objects }
    }

    fn world_bounds(&self) -> Bounds3 {
        unimplemented!()
    }

    fn intersect(&self, r: &Ray) -> Option<SurfaceInteraction> {
        let mut hit_anything: Option<SurfaceInteraction> = None;
        let mut closest_so_far = Float::INFINITY;

        for object in self.objects.iter() {
            if let Some(isect) = object.intersect(r) {
                if isect.dist < closest_so_far {
                    closest_so_far = isect.dist;
                    hit_anything = Some(isect);
                }
            }
        }
        return hit_anything;
    }

    fn intersect_p(&self, r: &Ray) -> bool {
        unimplemented!()
    }
}
