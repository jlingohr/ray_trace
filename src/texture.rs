use crate::interaction::SurfaceInteraction;

#[derive(Debug, Copy, Clone)]
pub enum Texture<T> {
    Constant(ConstantTexture<T>),
}

impl<T> Texture<T> {
    pub fn evaluate(&self, isect: &SurfaceInteraction) -> &T {
        match self {
            Texture::Constant(texture) => texture.evaluate(isect),
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub struct ConstantTexture<T> {
    value: T,
}

impl<T> ConstantTexture<T> {
    pub fn new(value: T) -> ConstantTexture<T> {
        ConstantTexture { value }
    }

    pub fn evaluate(&self, _isect: &SurfaceInteraction) -> &T {
        &self.value
    }
}
