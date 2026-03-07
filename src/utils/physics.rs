use crate::particle::particle::{Particle, ParticleType};
use crate::utils::operations::log_polynomial;

// Get particle energy
pub fn energy(particle: &Particle) -> f64 {
    (particle.state.p.mag().powf(2.0) + particle.state.m.powf(2.0)).sqrt()
}

// Get gamma factor of particle
pub fn gamma(particle: &Particle) -> f64 {
    let energy = energy(&particle);
    match particle.species {
        ParticleType::Gamma => panic!("gamma factor undefined for photons"),
        _                   => energy / particle.state.m,
    }
}

// Get beta factor of particle
pub fn beta(particle: &Particle) -> f64 {
    match particle.species {
        ParticleType::Gamma => 1.0,
        _                   => (1.0 - 1.0/gamma(particle).powf(2.0)).sqrt(),
    }
}

// Get dE/dx of ionizing particles (MeV/mm, hence the division by 10)
pub fn dEdx(particle: &Particle) -> f64 {
    let momentum = particle.state.p.mag();
    let value = match particle.species {
        ParticleType::Electron if momentum < 0.103 => 8.0 * 0.1,  // constant energy loss below fit range
        ParticleType::Electron                     => log_polynomial(momentum, vec![1.97185875, -4.90322067e-01, 5.67984147e-01, -3.78515229e-01, 1.96937857e-01, -6.69875048e-02, 1.30714285e-02, -1.31646064e-03, 5.29555090e-05]) * 0.1,
        ParticleType::Muon if momentum < 8.9       => 8.0 * 0.1,  // constant energy loss below fit range
        ParticleType::Muon if momentum < 50.0      => log_polynomial(momentum, vec![-2.21192313e+05, 4.16349323e+05, -3.02334049e+05, 9.22330794e+04, 1.78846389e+03, -9.81957228e+03, 2.97223872e+03, -3.90203242e+02, 1.99344973e+01]) * 0.1,
        ParticleType::Muon                         => log_polynomial(momentum, vec![1.13754387e+03, -1.13642381e+03, 4.96588219e+02, -1.23563655e+02, 1.91190645e+01, -1.88126582e+00, 1.14850292e-01, -3.97495919e-03, 5.96940644e-05]) * 0.1,
        ParticleType::Gamma                        => panic!("dE/dx is not defined for gammas"),
    };
    value
}

// Get kinetic energy of a particle in MeV
pub fn ke(particle: &Particle) -> f64 {
    let p = particle.state.p.mag();
    let m = particle.state.m;
    (p*p + m*m).sqrt() - m
}

// Tests
#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::vec3::Vec3;
    use approx::assert_relative_eq;

    fn assert_gamma_eq(lhs: Option<f64>, rhs: Option<f64>) {
        match (lhs, rhs) {
            (None, None)       => (),
            (None, Some(_))    => panic!("Mismatched variants: LHS is None, RHS is Some"),
            (Some(_), None)    => panic!("Mismatched variants: LHS is Some, RHS is None"),
            (Some(a), Some(b)) => assert_relative_eq!(a, b),
        }
    }

    #[test]
    fn test_physics_energy() {
        let p1 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0), ParticleType::Electron);
        let p2 = Particle::new(Vec3(6.6, -3.0, 5.2), Vec3(-3.0, 0.0, 4.0), ParticleType::Muon);
        let p3 = Particle::new(Vec3(-10.2, -3.5, -1.9), Vec3(10.0, -20.0, 30.0), ParticleType::Gamma);
        let p4 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0), ParticleType::Electron);
        let p5 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0), ParticleType::Muon);
        let p6 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0), ParticleType::Gamma);
        assert_relative_eq!(energy(&p1), 1.1229964381065507);
        assert_relative_eq!(energy(&p2), 105.77823783746825);
        assert_relative_eq!(energy(&p3), 37.416573867739416);
        assert_relative_eq!(energy(&p4), p4.state.m);
        assert_relative_eq!(energy(&p5), p5.state.m);
        assert_relative_eq!(energy(&p6), p6.state.m);
    }

    #[test]
    fn test_physics_gamma() {
        let p1 = Particle::new(Vec3(1.0, -1.0, 0.0), Vec3(4.0, -3.0, 0.0), ParticleType::Electron);
        let p2 = Particle::new(Vec3(0.0, 3.0, 5.0), Vec3(5.6, -2.1, -1.3), ParticleType::Muon);
        let p3 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0), ParticleType::Electron);
        let p4 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0), ParticleType::Muon);
        assert_relative_eq!(gamma(&p1), 9.835703071628355);
        assert_relative_eq!(gamma(&p2), 1.001676303733957);
        assert_relative_eq!(gamma(&p3), 1.0);
        assert_relative_eq!(gamma(&p4), 1.0);
        
    }

    #[test]
    fn test_physics_beta() {
        let p1 = Particle::new(Vec3(0.1, -1.5, -5.5), Vec3(-4.0, 0.0, 3.0), ParticleType::Electron);
        let p2 = Particle::new(Vec3(0.7, 9.8, 1.3), Vec3(3.0, 4.0, 0.0), ParticleType::Muon);
        let p3 = Particle::new(Vec3(1.5, -2.1, -4.8), Vec3(3.0, 4.0, 0.0), ParticleType::Gamma);
        assert_relative_eq!(beta(&p1), 0.9948181376436321);
        assert_relative_eq!(beta(&p2), 0.04726870197708133);
        assert_relative_eq!(beta(&p3), 1.0);
    }

    #[test]
    fn test_physics_dEdx() {
        let p1 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0), ParticleType::Electron);
        let p2 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0), ParticleType::Muon);
        let p3 = Particle::new(Vec3(1.0, -1.0, 0.0), Vec3(4.0, -3.0, 0.0), ParticleType::Electron);
        let p4 = Particle::new(Vec3(0.0, 3.0, 5.0), Vec3(5.6, -2.1, -1.3), ParticleType::Muon);
        let p5 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(9.7, 15.2, 51.1), ParticleType::Electron);
        let p6 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(53.4, -98.3, -89.5), ParticleType::Muon);
        assert_relative_eq!(dEdx(&p1), 0.8);
        assert_relative_eq!(dEdx(&p2), 0.8);
        assert_relative_eq!(dEdx(&p3), 0.18667002945559819);
        assert_relative_eq!(dEdx(&p4), 0.8);
        assert_relative_eq!(dEdx(&p5), 0.21359254760465154);
        assert_relative_eq!(dEdx(&p6), 0.24890235819417583);
    }

    #[test]
    fn test_physics_ke() {
        let p1 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0), ParticleType::Electron);
        let p2 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0), ParticleType::Muon);
        let p3 = Particle::new(Vec3(1.0, -1.0, 0.0), Vec3(4.0, -3.0, 0.0), ParticleType::Electron);
        let p4 = Particle::new(Vec3(0.0, 3.0, 5.0), Vec3(4.0, -3.0, 0.0), ParticleType::Muon);
        let p5 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(9.7, 15.2, 51.1), ParticleType::Electron);
        let p6 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(53.4, -98.3, -89.5), ParticleType::Muon);
        assert_relative_eq!(ke(&p1), 0.0);
        assert_relative_eq!(ke(&p2), 0.0);
        assert_relative_eq!(ke(&p3), 4.51504426960209);
        assert_relative_eq!(ke(&p4), 0.11823783746825711);
        assert_relative_eq!(ke(&p5), 53.679415397928075);
        assert_relative_eq!(ke(&p6), 72.35330175017819);
    }
}
