use crate::utils::vec3::Vec3;

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

// Will have a `track` field, and probably a `interaction_dist` field
pub struct Photon {
    pub state: PhotonState,
}

impl Photon {
    pub fn new(position: Vec3, direction: Vec3, time: f64) -> Self {
        let photon_state = PhotonState::new(position, direction, time);
        Photon { state: photon_state }
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
}
