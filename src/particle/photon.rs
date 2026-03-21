use rand::Rng;
use rand_distr::{Distribution, UnitSphere};

use crate::utils::vec3::Vec3;
use crate::geometry::volume::Volume;
use crate::utils::constants::C;

pub struct PhotonState {
    pub r: Vec3,  // position (mm)
    pub d: Vec3,  // direction (1)
    pub t: f64,   // time (ns)
    pub alive: bool,
}

impl PhotonState {
    pub fn new(position: Vec3, direction: Vec3, time: f64) -> Self {
        PhotonState { r: position, d: direction, t: time, alive: true }
    }
}

pub struct PhotonTrack {
    pub vertices: Vec<(Vec3, f64)>,  // (position, time) of each vertex
}

impl PhotonTrack {
    pub fn new(origin: Vec3, time: f64) -> Self {
        PhotonTrack { vertices: vec![(origin, time)] }
    }

    pub fn record(&mut self, position: Vec3, time: f64) {
        self.vertices.push((position, time));
    }
}

// Will probably have a `interaction_dist` field
pub struct Photon {
    pub state: PhotonState,
    pub track: PhotonTrack,
}

impl Photon {
    pub fn new(position: Vec3, direction: Vec3, time: f64) -> Self {
        let photon_state = PhotonState::new(position, direction, time);
        let photon_track = PhotonTrack::new(position, time);
        Photon { state: photon_state, track: photon_track }
    }

    pub fn simulate(&mut self, volume: &Volume, rng: &mut impl Rng) {
        let photon_speed = C / volume.n;

        loop {
            // Check the three possible interaction distances
            let detector_dist = volume.intersect(self.state.r, self.state.d);
            let absorption_dist = -volume.la * rng.random::<f64>().ln();
            let scattering_dist = -volume.ls * rng.random::<f64>().ln();

            // Record the new photon vertex
            let distance = f64::min(f64::min(detector_dist, absorption_dist), scattering_dist);
            self.state.r += self.state.d * distance;
            self.state.t += distance / photon_speed;
            self.record(self.state.r, self.state.t);
            
            // Perform correct behaviour depending on interaction type
            if scattering_dist < detector_dist && scattering_dist < absorption_dist { // scatter
                let [x, y, z] = UnitSphere.sample(rng);
                self.state.d = Vec3(x, y, z);
            } else { // detector wall collision or absorption
                break;
            }
        }
    }

    pub fn record(&mut self, position: Vec3, time: f64) {
        self.track.record(position, time);
    }
}


// Tests
#[cfg(test)]
mod tests {
    use super::*;
    use crate::assert_vec3_eq;
    use approx::assert_relative_eq;

    #[test]
    fn test_photon_creation() {
        let v1 = Vec3(1.2, 4.3, -2.2);
        let v2 = Vec3(0.8, -3.3, 7.1);
        let d1 = Vec3(-3.5, 2.1, -5.0);
        let d2 = Vec3(4.4, 2.9, -15.0);
        let t1 = 0.5;
        let t2 = 0.97;
        let ph1 = Photon::new(v1, d1, t1);
        let ph2 = Photon::new(v2, d2, t2);
        assert_vec3_eq!(ph1.state.r, v1);
        assert_vec3_eq!(ph2.state.r, v2);
        assert_vec3_eq!(ph1.state.d, d1);
        assert_vec3_eq!(ph2.state.d, d2);
        assert_relative_eq!(ph1.state.t, t1);
        assert_relative_eq!(ph2.state.t, t2);
    }

    #[test]
    fn test_photon_record() {
        let v1 = Vec3(1.2, 4.3, -2.2);
        let v2 = Vec3(0.8, -3.3, 7.1);
        let v3 = Vec3(-3.5, 2.1, -5.0);
        let v4 = Vec3(4.4, 2.9, -15.0);
        let d = Vec3(1.0, 0.0, 0.0);
        let mut ph = Photon::new(v1, d, 0.0);
        assert_vec3_eq!(ph.track.vertices[0].0, v1);
        assert_relative_eq!(ph.track.vertices[0].1, 0.0);
        ph.record(v2, 0.5);
        assert_vec3_eq!(ph.track.vertices[1].0, v2);
        assert_relative_eq!(ph.track.vertices[1].1, 0.5);
        ph.record(v3, 0.9);
        assert_vec3_eq!(ph.track.vertices[2].0, v3);
        assert_relative_eq!(ph.track.vertices[2].1, 0.9);
        ph.record(v4, 1.7);
        assert_vec3_eq!(ph.track.vertices[3].0, v4);
        assert_relative_eq!(ph.track.vertices[3].1, 1.7);
    }
}
