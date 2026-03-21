use crate::particle::particle::Particle;
use crate::utils::vec3::Vec3;

#[allow(non_snake_case)]
#[derive(Clone)]
pub struct Volume {
    pub size: f64,  // cube edge length (mm)
    pub X0: f64,    // radiation length (mm)
    pub LY: f64,    // light yield (photons / MeV)
    pub n: f64,     // refractive index
    pub la: f64,    // absorption length (mm)
    pub ls: f64,    // scattering length (mm)
}

impl Volume {
    pub fn new(s: f64, radiation_length: f64, light_yield: f64, refractive_index: f64, absorption_length: f64, scattering_length: f64) -> Self {
        Volume { size: s, X0: radiation_length, LY: light_yield, n: refractive_index, la: absorption_length, ls: scattering_length }
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

    // Calculate the intersection distance of a ray with the volume
    // Returns the distance from the ray origin to the intersection point
    // Assumes a normalized direction
    pub fn intersect(&self, position: Vec3, direction: Vec3) -> f64 {
        let hs = self.size / 2.0;  // half-size
        let mut t_min = f64::INFINITY;

        // Solve:  t_i = (+-hs - r_i) / d_i  for each axis
        // The smallest t is selected

        // Check each axis
        for i in 0..3 {
            let r = [position.0, position.1, position.2][i];
            let d = [direction.0, direction.1, direction.2][i];

            // Skip if ray is parallel to this pair of walls
            if d.abs() < 1e-10 {
                continue;
            }

            // Check both walls along this axis
            for &wall in &[hs, -hs] {
                // Solve for t
                let t = (wall - r) / d;
                // Only consider intersections in front of the ray
                if t > 0.0 && t < t_min {
                    t_min = t;
                }
            }
        }

        t_min
    }
}


// Tests
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use crate::assert_vec3_eq;
    use crate::particle::particle::ParticleType;

    #[test]
    fn test_volume_creation() {
        let v1 = Volume::new(5.0, 50.0, 200.0, 1.2, 2000.0, 2000.0);
        let v2 = Volume::new(15.0, 36.0, 150.0, 1.5, 2500.0, 1500.0);
        let v3 = Volume::new(62.3, 52.0, 300.0, 1.4, 1500.0, 1000.0);
        assert_relative_eq!(v1.size, 5.0);
        assert_relative_eq!(v2.size, 15.0);
        assert_relative_eq!(v3.size, 62.3);
    }

    #[test]
    fn test_volume_contains() {
        let v1 = Volume::new(10.0, 30.0, 110.0, 1.7, 1000.0, 800.0);
        let v2 = Volume::new(28.4, 44.2, 130.0, 1.0, 1200.0, 1500.0);
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

    #[test]
    fn test_intersect() {
        let v1 = Volume::new(100.0, 30.0, 110.0, 1.3, 2000.0, 1800.0);
        let v2 = Volume::new(150.0, 44.2, 130.0, 1.54, 5000.0, 3000.0);
        let p1 = Vec3::new(0.0, 0.0, 0.0);
        let p2 = Vec3::new(40.0, 30.0, -20.0);
        let d1 = Vec3::new(1.0, 0.0, 0.0);
        let d2 = Vec3::new(0.0, -1.0, 0.0);
        let d3 = Vec3::new(0.2672612419124243, 0.5345224838248487, 0.8017837257372731);
        assert_relative_eq!(v1.intersect(p1, d1), 50.0);
        assert_relative_eq!(v1.intersect(p2, d1), 10.0);
        assert_relative_eq!(v1.intersect(p1, d2), 50.0);
        assert_relative_eq!(v1.intersect(p2, d2), 80.0);
        assert_relative_eq!(v1.intersect(p1, d3), 62.360956446232365);
        assert_relative_eq!(v1.intersect(p2, d3), 37.41657386773942);
        assert_relative_eq!(v2.intersect(p1, d1), 75.0);
        assert_relative_eq!(v2.intersect(p2, d1), 35.0);
        assert_relative_eq!(v2.intersect(p1, d2), 75.0);
        assert_relative_eq!(v2.intersect(p2, d2), 105.0);
        assert_relative_eq!(v2.intersect(p1, d3), 93.54143466934855);
        assert_relative_eq!(v2.intersect(p2, d3), 84.18729120241369);
    }
}
