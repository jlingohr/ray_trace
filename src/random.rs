use hexf::*;
use rand::rngs::ThreadRng;
use rand::Rng;

use crate::pbrt::Float;

pub const FLOAT_ONE_MINUS_EPSILON: Float = hexf64!("0x1.fffffep-1");
//TODO this is for f32
//#endif
pub const PCG32_DEFAULT_STATE: u64 = 0x853c_49e6_748f_ea9b;
pub const PCG32_DEFAULT_STREAM: u64 = 0xda3e_39cb_94b9_5bdb;
pub const PCG32_MULT: u64 = 0x5851_f42d_4c95_7f2d;

pub struct Random {
    state: u64,
    inc: u64,
}

impl Random {
    pub fn new() -> Random {
        Random {
            state: PCG32_DEFAULT_STATE,
            inc: PCG32_DEFAULT_STREAM,
        }
    }

    pub fn set_sequence(&mut self, initseq: u64) {
        self.state = 0_u64;
        let (shl, _overflow) = initseq.overflowing_shl(1);
        self.inc = shl | 1;
        self.uniform_uint32();
        let (add, _overflow) = self.state.overflowing_add(PCG32_DEFAULT_STATE);
        self.state = add;
        self.uniform_uint32();
    }
    pub fn uniform_uint32(&mut self) -> u32 {
        let oldstate: u64 = self.state;
        // C++: state = oldstate * PCG32_MULT + inc;
        let (mul, _overflow) = oldstate.overflowing_mul(PCG32_MULT);
        let (add, _overflow) = mul.overflowing_add(self.inc);
        self.state = add;
        // C++: uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
        let (shr, _overflow) = oldstate.overflowing_shr(18);
        let combine = shr ^ oldstate;
        let (shr, _overflow) = combine.overflowing_shr(27);
        let xorshifted: u32 = shr as u32;
        // C++: uint32_t rot = (uint32_t)(oldstate >> 59u);
        let (shr, _overflow) = oldstate.overflowing_shr(59);
        let rot: u32 = shr as u32;
        // C++: return (xorshifted >> rot) | (xorshifted << ((~rot + 1u) & 31));
        let (shr, _overflow) = xorshifted.overflowing_shr(rot);
        // bitwise not in Rust is ! (not the ~ operator like in C)
        let neg = !rot;
        let (add, _overflow) = neg.overflowing_add(1_u32);
        let (shl, _overflow) = xorshifted.overflowing_shl(add & 31);
        shr | shl
    }
    pub fn uniform_uint32_bounded(&mut self, b: u32) -> u32 {
        // bitwise not in Rust is ! (not the ~ operator like in C)
        let threshold = (!b + 1) & b;
        loop {
            let r = self.uniform_uint32();
            if r >= threshold {
                return r % b;
            }
        }
    }
    pub fn uniform_float(&mut self) -> Float {
        (self.uniform_uint32() as Float * hexf32!("0x1.0p-32") as Float)
            .min(FLOAT_ONE_MINUS_EPSILON)
    }
}
