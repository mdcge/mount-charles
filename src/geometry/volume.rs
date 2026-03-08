use crate::particle::particle::Particle;
use crate::utils::vec3::Vec3;

#[allow(non_snake_case)]
#[derive(Clone)]
pub struct Volume {
    pub size: f64,  // cube edge length (mm)
    pub X0: f64,    // radiation length (mm)
}

impl Volume {
    pub fn new(s: f64, rad_len: f64) -> Self {
        Volume { size: s, X0: rad_len }
    }

    pub fn contains(&self, particle: &Particle) -> bool {
        let Vec3(x, y, z) = particle.state.r;
        let hs = self.size / 2.0;
        if (-hs <= x && x <= hs) && (-hs <= y && y <= hs) && (-hs <= z && z <= hs) {
            true
        } else {
            false
        }
    }
}


// Tests
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use crate::particle::particle::ParticleType;

    #[test]
    fn test_volume_creation() {
        let v1 = Volume::new(5.0, 50.0);
        let v2 = Volume::new(15.0, 36.0);
        let v3 = Volume::new(62.3, 52.0);
        assert_relative_eq!(v1.size, 5.0);
        assert_relative_eq!(v2.size, 15.0);
        assert_relative_eq!(v3.size, 62.3);
    }

    #[test]
    fn test_volume_contains() {
        let v1 = Volume::new(10.0, 30.0);
        let v2 = Volume::new(28.4, 44.2);
        let p1 = Particle::new(Vec3(1.0, 2.0, -3.0), Vec3(5.0, 0.0, 0.0), ParticleType::Electron);
        let p2 = Particle::new(Vec3(4.2, -1.5, 5.1), Vec3(3.4, -2.0, 0.7), ParticleType::Muon);
        let p3 = Particle::new(Vec3(20.1, -10.3, -9.7), Vec3(-100.0, 0.0, -52.1), ParticleType::Gamma);
        assert_eq!(v1.contains(&p1), true);
        assert_eq!(v1.contains(&p2), false);
        assert_eq!(v1.contains(&p3), false);
        assert_eq!(v2.contains(&p1), true);
        assert_eq!(v2.contains(&p2), true);
        assert_eq!(v2.contains(&p3), false);
    }
}
