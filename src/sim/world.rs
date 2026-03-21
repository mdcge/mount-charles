use rand::{Rng, rngs::StdRng, SeedableRng};
use rand_distr::{Distribution, UnitSphere};
use rayon::prelude::*;

use crate::particle::particle::{Particle, ParticleType};
use crate::particle::photon::{Photon, PhotonTrack};
use crate::particle::track::Track;
use crate::geometry::volume::Volume;
use crate::utils::physics::{ke, lambda};
use crate::utils::vec3::Vec3;

pub struct World {
    pub time: f64,    // world time (ns)
    pub dt: f64,  // time step (ns)
    pub particles: Vec<Particle>,
    pub photons: Vec<Photon>,
    pub volume: Volume,
    pub rng: StdRng,
}

impl World {
    pub fn new(particle_list: Vec<Particle>, vol: Volume, timestep: f64, random_seed: u64) -> Self {
        World {
            time: 0.0,
            dt: timestep,
            particles: particle_list,
            photons: vec![],
            volume: vol,
            rng: StdRng::seed_from_u64(random_seed)
        }
    }

    pub fn has_alive_particles(&self) -> bool {
        self.particles.iter().any(|p| p.state.alive)
    }

    pub fn step(&mut self) {
        // Increment time
        self.time += self.dt;
        
        for particle in &mut self.particles {            
            // Check if particle KE is below 10keV (if not a gamma)
            if ke(&particle) < 0.01 && !matches!(particle.species, ParticleType::Gamma) {
                particle.state.alive = false;
                continue;
            }

            // Ignore if particle is dead
            if !particle.state.alive {
                continue;
            }

            // Propagate the particle
            particle.propagate(self.dt);

            // Check if particle is out of bounds
            if !self.volume.contains(&particle) {
                particle.state.alive = false;
                continue;
            }

            // Energy before interaction
            let pre_energy = ke(&particle);

            // Interact the particle if distance to interaction is 0
            if particle.interaction_dist <= 0.0 {
                particle.interact(&mut self.rng, self.volume.X0, self.dt);
                particle.interaction_dist = match particle.species {
                    ParticleType::Electron | ParticleType::Muon => 0.0,  // forces interaction every step
                    ParticleType::Gamma                         => -lambda(particle) * f64::ln(self.rng.random::<f64>()),
                }
            }

            // Energy after interaction
            let post_energy = ke(&particle);
            let energy_deposit = match pre_energy - post_energy {
                0.0 => None,
                #[allow(non_snake_case)]
                dE  => Some(f64::min(dE, pre_energy)),
            };

            // Record track point
            particle.track.record(particle.state.r, self.time, energy_deposit);
        }
    }

    pub fn simulate_scintillation(&mut self) {
        // Collect energy deposits from all particles
        let energy_deposits = self.particles
                                  .iter()
                                  .flat_map(|p| p.track.points.iter())
                                  .filter_map(|point| {
                                      point.E.map(|energy| (point.r, point.t, energy))
                                  })
                                  .collect::<Vec<(Vec3, f64, f64)>>();
        
        // Generate photons for each energy deposit
        for (position, time, energy) in energy_deposits {
            let nb_photons = (energy * self.volume.LY) as u64;
            for _ in 0..nb_photons {
                let [x, y, z] = UnitSphere.sample(&mut self.rng);
                let direction = Vec3(x, y, z);
                
                let mut photon = Photon::new(position, direction, time);
                photon.simulate(&self.volume);
                
                self.photons.push(photon);
            }
        }
    }

    pub fn tracks(&self) -> Vec<&Track> {
        self.particles.iter().map(|p| &p.track).collect()
    }

    pub fn photon_tracks(&self) -> Vec<&PhotonTrack> {
        self.photons.iter().map(|ph| &ph.track).collect()
    }
}


// Tests
#[cfg(test)]
mod tests{
    use super::*;
    use crate::particle::particle::ParticleType;
    use crate::utils::vec3::Vec3;

    #[test]
    fn test_world_creation() {
        let v1 = Volume::new(10.0, 53.2, 200.0);
        let v2 = Volume::new(28.4, 60.0, 250.0);
        let p1 = Particle::new(Vec3(1.0, 2.0, -3.0), Vec3(5.0, 0.0, 0.0), ParticleType::Electron);
        let p2 = Particle::new(Vec3(4.2, -1.5, 5.1), Vec3(3.4, -2.0, 0.7), ParticleType::Muon);
        let p3 = Particle::new(Vec3(20.1, -10.3, -9.7), Vec3(-100.0, 0.0, -52.1), ParticleType::Gamma);
        let w1 = World::new(vec![p1, p2, p3], v1.clone(), 0.1, 0);

        let p1 = Particle::new(Vec3(1.0, 2.0, -3.0), Vec3(5.0, 0.0, 0.0), ParticleType::Electron);
        let p3 = Particle::new(Vec3(20.1, -10.3, -9.7), Vec3(-100.0, 0.0, -52.1), ParticleType::Gamma);
        let w2 = World::new(vec![p1, p3], v2.clone(), 0.01, 1);

        let p2 = Particle::new(Vec3(4.2, -1.5, 5.1), Vec3(3.4, -2.0, 0.7), ParticleType::Muon);
        let w3 = World::new(vec![p2], v2, 0.005, 2);
        let w4 = World::new(vec![], v1, 1.0, 3);
        assert_eq!(w1.particles.len(), 3);
        assert_eq!(w2.particles.len(), 2);
        assert_eq!(w3.particles.len(), 1);
        assert_eq!(w4.particles.len(), 0);
    }

    #[test]
    fn test_world_has_alive_particles() {
        let v1 = Volume::new(10.0, 53.2, 135.0);
        let v2 = Volume::new(28.4, 60.0, 96.0);
        let p1 = Particle::new(Vec3(1.0, 2.0, -3.0), Vec3(5.0, 0.0, 0.0), ParticleType::Electron);
        let p2 = Particle::new(Vec3(4.2, -1.5, 5.1), Vec3(3.4, -2.0, 0.7), ParticleType::Muon);
        let p3 = Particle::new(Vec3(20.1, -10.3, -9.7), Vec3(-100.0, 0.0, -52.1), ParticleType::Gamma);
        let w1 = World::new(vec![p1, p2, p3], v1.clone(), 0.1, 15);
        let w2 = World::new(vec![], v2.clone(), 0.01, 837);

        let p2 = Particle::new(Vec3(4.2, -1.5, 5.1), Vec3(3.4, -2.0, 0.7), ParticleType::Muon);
        let w3 = World::new(vec![p2], v2, 0.005, 9882);
        let w4 = World::new(vec![], v1, 1.0, 21);
        assert!(w1.has_alive_particles());
        assert!(!w2.has_alive_particles());
        assert!(w3.has_alive_particles());
        assert!(!w4.has_alive_particles());
    }
}
