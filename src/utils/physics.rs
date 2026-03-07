use crate::particle::particle::{Particle, ParticleType};
use crate::utils::operations::log_polynomial;
use crate::utils::constants::{ELECTRON_DEDX_COEFFS, MUON_DEDX_LOW_COEFFS, MUON_DEDX_HIGH_COEFFS, COMPTON_COEFFS, PHOTOELECTRIC_A_COEFFS, PHOTOELECTRIC_P_COEFFS};

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
        ParticleType::Electron                     => log_polynomial(momentum, ELECTRON_DEDX_COEFFS.into()) * 0.1,
        ParticleType::Muon if momentum < 8.9       => 8.0 * 0.1,  // constant energy loss below fit range
        ParticleType::Muon if momentum < 50.0      => log_polynomial(momentum, MUON_DEDX_LOW_COEFFS.into()) * 0.1,
        ParticleType::Muon                         => log_polynomial(momentum, MUON_DEDX_HIGH_COEFFS.into()) * 0.1,
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

// Compton scattering attenuation coefficient
fn mu_compton(particle: &Particle) -> f64 {
    log_polynomial(particle.state.p.mag(), COMPTON_COEFFS.into())
}

// Photoelectric effect attenuation coefficient
fn mu_photo(particle: &Particle) -> f64 {
    let energy = particle.state.p.mag();
    PHOTOELECTRIC_A_COEFFS[0] * energy.powf(-PHOTOELECTRIC_P_COEFFS[0]) + PHOTOELECTRIC_A_COEFFS[1] * energy.powf(-PHOTOELECTRIC_P_COEFFS[1])
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

    #[test]
    fn test_mu_compton() {
        let p1 = Particle::new(Vec3(1.5, -2.1, -4.8), Vec3(3.0, 4.0, 0.0), ParticleType::Gamma);
        let p2 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(-10.0, -100.0, 20.0), ParticleType::Gamma);
        let p3 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(0.5, 0.0, 0.0), ParticleType::Gamma);
        assert_relative_eq!(mu_compton(&p1), 0.027507314428759033);
        assert_relative_eq!(mu_compton(&p2), 0.002409009744319685);
        assert_relative_eq!(mu_compton(&p3), 0.09734846035725048);
    }

    #[test]
    fn test_mu_photo() {
        let p1 = Particle::new(Vec3(1.5, -2.1, -4.8), Vec3(3.0, 4.0, 0.0), ParticleType::Gamma);
        let p2 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(-10.0, -100.0, 20.0), ParticleType::Gamma);
        let p3 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(0.5, 0.0, 0.0), ParticleType::Gamma);
        assert_relative_eq!(mu_photo(&p1), 2.9890205789156227e-7);
        assert_relative_eq!(mu_photo(&p2), 1.2436722113325266e-8);
        assert_relative_eq!(mu_photo(&p3), 2.1008095955227985e-5);
    }
}
