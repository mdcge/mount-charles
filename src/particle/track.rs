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
    }
}
