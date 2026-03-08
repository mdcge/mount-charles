use crate::utils::vec3::Vec3;
use crate::particle::particle::ParticleType;

// Single track point
pub struct TrackPoint {
    pub r: Vec3,         // position
    pub t: f64,          // time
    pub E: Option<f64>,  // energy deposit (`None` if no energy deposited)
}

impl TrackPoint {
    pub fn new(position: Vec3, time: f64, energy_deposit: Option<f64>) -> Self {
        TrackPoint { r: position, t: time, E: energy_deposit }
    }
}


// Full track
pub struct Track {
    pub particle_type: ParticleType,
    pub points: Vec<TrackPoint>,
}

impl Track {
    pub fn new(particle_type: ParticleType, position: Vec3, time: f64) -> Self {
        Track { particle_type: particle_type, points: vec![TrackPoint::new(position, time, None)] }
    }

    pub fn new_empty(particle_type: ParticleType) -> Self {
        Track { particle_type: particle_type, points: vec![] }
    }

    pub fn record(&mut self, position: Vec3, time: f64, energy_deposit: Option<f64>) {
        self.points.push(TrackPoint::new(position, time, energy_deposit));
    }
}


// Tests
#[cfg(test)]
mod tests {
    use super::*;
    use crate::assert_vec3_eq;
    use approx::assert_relative_eq;

    #[test]
    fn test_trackpoint_creation() {
        let tp1 = TrackPoint::new(Vec3(0.0, 0.0, 0.0), 0.0, None);
        let tp2 = TrackPoint::new(Vec3(1.7, -3.2, -5.0), 4.2, Some(13.5));
        assert_vec3_eq!(tp1.r, Vec3(0.0, 0.0, 0.0));
        assert_relative_eq!(tp1.t, 0.0);
        assert_eq!(tp1.E, None);
        assert_vec3_eq!(tp2.r, Vec3(1.7, -3.2, -5.0));
        assert_relative_eq!(tp2.t, 4.2);
        assert_eq!(tp2.E, Some(13.5));
    }

    #[test]
    fn test_track_creation() {
        let _t1 = Track::new(ParticleType::Electron, Vec3(0.0, 0.0, 0.0), 0.0);
        let _t2 = Track::new(ParticleType::Gamma, Vec3(-7.3, 5.2, -10.1), 5.9);
    }

    #[test]
    fn test_track_record() {
        let mut t1 = Track::new(ParticleType::Muon, Vec3(0.0, 0.0, 0.0), 0.0);
        t1.record(Vec3(1.0, 0.0, 0.0), 1.0, None);
        t1.record(Vec3(2.5, -0.3, 0.0), 2.0, Some(5.0));
        t1.record(Vec3(3.1, -0.8, 0.0), 2.5, None);
        assert_vec3_eq!(t1.points[0].r, Vec3(0.0, 0.0, 0.0));
        assert_relative_eq!(t1.points[0].t, 0.0);
        assert_eq!(t1.points[0].E, None);
        assert_vec3_eq!(t1.points[1].r, Vec3(1.0, 0.0, 0.0));
        assert_relative_eq!(t1.points[1].t, 1.0);
        assert_eq!(t1.points[1].E, None);
        assert_vec3_eq!(t1.points[2].r, Vec3(2.5, -0.3, 0.0));
        assert_relative_eq!(t1.points[2].t, 2.0);
        assert_eq!(t1.points[2].E, Some(5.0));
        assert_vec3_eq!(t1.points[3].r, Vec3(3.1, -0.8, 0.0));
        assert_relative_eq!(t1.points[3].t, 2.5);
        assert_eq!(t1.points[3].E, None);
    }
}
