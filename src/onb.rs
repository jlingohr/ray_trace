use super::vector::Vec3;

pub struct ONB {
    axis: [Vec3; 3],
}

impl ONB {
    pub fn build_from_w(n: Vec3) -> ONB {
        let w = n.unit();
        let a = if w[0].abs() > 0.9 {
            Vec3::new(0.0, 1.0, 0.0)
        } else {
            Vec3::new(1.0, 0.0, 0.0)
        };
        let v = w.cross(&a).unit();
        let u = w.cross(&v);

        let axis = [u, v, w];
        ONB { axis }
    }

    pub fn u(&self) -> Vec3 {
        self.axis[0]
    }

    pub fn v(&self) -> Vec3 {
        self.axis[1]
    }

    pub fn w(&self) -> Vec3 {
        self.axis[2]
    }

    pub fn local(&self, a: Vec3) -> Vec3 {
        a.x * self.u() + a.y * self.v() + a.z * self.w()
    }
}
