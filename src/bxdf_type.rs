use bitflags::bitflags;

// For each BxDF the flags should have at least one of REFLECTION or TRANSMISSION set
// and exactly one of diffuse, glossy, and specular flags.
bitflags! {
    pub struct BxDFType: u8 {
        const Empty = 0b00000000;
        const BsdfReflection = 0b00000001;
        const BsdfTransmission = 0b00000010;
        const BsdfDiffuse = 0b00000100;
        const BsdfGlossy = 0b00001000;
        const BsdfSpecular = 0b00010000;
        const BsdfAll = 0b11111111;
    }
}
