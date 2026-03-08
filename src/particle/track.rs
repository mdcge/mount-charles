use crate::utils::vec3::Vec3;

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
    pub points: Vec<TrackPoint>,
}

impl Track {
    pub fn new(position: Vec3, time: f64) -> Self {
        Track { points: vec![TrackPoint::new(position, time, None)] }
    }

    pub fn new_empty() -> Self {
        Track { points: vec![] }
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
}
