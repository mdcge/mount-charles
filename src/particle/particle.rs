use rand::Rng;

use crate::utils::vec3::Vec3;
use crate::utils::constants::{C, Me, Mmu, Mg};
use crate::utils::physics::{beta, dEdx, ke, mu_compton, mu_total, compton_angle};
use crate::particle::track::Track;

// Particle state
#[derive(Debug, Clone)]
pub struct ParticleState {
    pub r: Vec3,  // position (mm)
    pub p: Vec3,  // momentum (MeV)
    pub m: f64,   // mass (MeV)
    pub alive: bool,
}

impl ParticleState {
    pub fn new(pos: Vec3, mom: Vec3, mass: f64) -> Self {
        ParticleState { r: pos, p: mom, m: mass, alive: true }
    }
}

// Particle type
#[derive(Clone)]
pub enum ParticleType {
    Electron,
    Muon,
    Gamma,
}

// Particle
pub struct Particle {
    pub species: ParticleType,
    pub state: ParticleState,
    pub track: Track,
    pub interaction_dist: f64,  // distance to next interaction
}

impl Particle {
    pub fn new(pos: Vec3, mom: Vec3, part_type: ParticleType) -> Self {
        let mass = match part_type {
            ParticleType::Electron => Me,
            ParticleType::Muon     => Mmu,
            ParticleType::Gamma    => Mg,
        };
        let particle_state = ParticleState::new(pos, mom, mass);

        Particle { species: part_type.clone(), state: particle_state, track: Track::new(part_type, pos, 0.0), interaction_dist: 0.0 }
    }

    pub fn propagate(&mut self, dt: f64) {
        let dir = self.state.p.norm();
        let beta = beta(&self);
        let dx = beta * C * dt;  // beta[1] * C[mm/ns] * dt[ns]
        self.state.r += dir * dx;  // dir[1, 1, 1] * dx[mm]
        self.interaction_dist -= dx;
    }

    #[allow(non_snake_case)]
    pub fn interact(&mut self, rng: &mut impl Rng, X0: f64, dt: f64) {
        match self.species {
            ParticleType::Electron | ParticleType::Muon => self.interact_ionising(rng, X0, dt),
            ParticleType::Gamma                         => self.interact_gamma(rng),
        }
    }

    #[allow(non_snake_case)]
    fn interact_ionising(&mut self, rng: &mut impl Rng, X0: f64, dt: f64) {
        let p = self.state.p.mag();
        let ke_pre = ke(&self);  // kinetic energy before step

        let beta = beta(&self);
        let dx = beta * C * dt;  // step size (mm)
        let dEdx = dEdx(&self);
        let theta0 = (13.6 / (beta * p)) * (dx/X0).sqrt() * (1.0 + 0.038 * f64::ln(dx/X0));

        // Subtract energy lost in step, clamp to 0 if negative
        let ke_post = f64::max(ke_pre - (dEdx * dx), 0.0);
        // Resize momentum vector with new momentum
        self.state.p = self.state.p.norm() * (ke_post * (ke_post + 2.0*self.state.m)).sqrt();
        // Deflect momentum vector
        self.state.p.deflect(rng, theta0);
    }

    fn interact_gamma(&mut self, rng: &mut impl Rng) {
        let mu_compton = mu_compton(&self);
        let mu_total = mu_total(&self);
        let p_compton = mu_compton / mu_total;
        if rng.random::<f64>() < p_compton {
            // Compton scattering
            let theta = compton_angle(&self, rng);
            let energy = self.state.p.mag();
            let new_energy = energy / (1.0 + (energy / Me) * (1.0 - f64::cos(theta)));
            self.state.p = self.state.p.norm() * new_energy;
            self.state.p.rotate(rng, theta);
        } else {
            // Photoelectric effect
            self.state.p = Vec3(0.0, 0.0, 0.0);
            self.state.alive = false;
        }
    }
}


// Tests
#[cfg(test)]
mod tests {
    use super::*;
    use crate::assert_vec3_eq;
    use approx::assert_relative_eq;

    #[test]
    fn test_particlestate_creation() {
        let v1 = Vec3(1.2, 4.3, -2.2);
        let v2 = Vec3(0.8, -3.3, 7.1);
        let m1 = 50.0;
        let m2 = 73.1;
        let ps1 = ParticleState::new(v1, v2, m1);
        let ps2 = ParticleState::new(v2, v1, m2);
        assert_vec3_eq!(ps1.r, v1);
        assert_vec3_eq!(ps1.p, v2);
        assert_relative_eq!(ps1.m, m1);
        assert_eq!(ps1.alive, true);
        assert_vec3_eq!(ps2.r, v2);
        assert_vec3_eq!(ps2.p, v1);
        assert_relative_eq!(ps2.m, m2);
        assert_eq!(ps2.alive, true);
    }

    #[test]
    fn test_particletype_creation() {
        let _pt1 = ParticleType::Electron;
        let _pt2 = ParticleType::Muon;
        let _pt3 = ParticleType::Gamma;
    }

    #[test]
    fn test_particle_creation() {
        let electron = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0), ParticleType::Electron);
        let muon = Particle::new(Vec3(5.0, -2.0, 10.0), Vec3(2.0, 3.0, -4.0), ParticleType::Muon);
        let gamma = Particle::new(Vec3(-1.2, 7.6, 6.7), Vec3(-9.8, -2.5, -1.1), ParticleType::Gamma);
        assert_vec3_eq!(electron.state.r, Vec3(0.0, 0.0, 0.0));
        assert_vec3_eq!(electron.state.p, Vec3(1.0, 0.0, 0.0));
        assert_relative_eq!(electron.state.m, 0.511);
        assert_eq!(electron.state.alive, true);
        assert_vec3_eq!(muon.state.r, Vec3(5.0, -2.0, 10.0));
        assert_vec3_eq!(muon.state.p, Vec3(2.0, 3.0, -4.0));
        assert_relative_eq!(muon.state.m, 105.66);
        assert_eq!(muon.state.alive, true);
        assert_vec3_eq!(gamma.state.r, Vec3(-1.2, 7.6, 6.7));
        assert_vec3_eq!(gamma.state.p, Vec3(-9.8, -2.5, -1.1));
        assert_relative_eq!(gamma.state.m, 0.0);
        assert_eq!(gamma.state.alive, true);
    }

    #[test]
    fn test_particle_propagate() {
        let mut electron1 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0), ParticleType::Electron);
        let mut muon1 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0), ParticleType::Muon);
        let mut gamma1 = Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 1.0), ParticleType::Gamma);

        electron1.propagate(1.0);
        assert_vec3_eq!(electron1.state.r, Vec3(266.9576214377587, 0.0, 0.0));
        electron1.propagate(0.1);
        assert_vec3_eq!(electron1.state.r, Vec3(266.9576214377587+26.69576214377587, 0.0, 0.0));
        muon1.propagate(1.0);
        assert_vec3_eq!(muon1.state.r, Vec3(0.0, 2.837204544727953, 0.0));
        muon1.propagate(0.1);
        assert_vec3_eq!(muon1.state.r, Vec3(0.0, 2.837204544727953+0.2837204544727953, 0.0));
        gamma1.propagate(1.0);
        assert_vec3_eq!(gamma1.state.r, Vec3(0.0, 0.0, 299.792458));
        gamma1.propagate(0.1);
        assert_vec3_eq!(gamma1.state.r, Vec3(0.0, 0.0, 299.792458+29.9792458));
    }
}
