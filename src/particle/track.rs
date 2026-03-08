use crate::utils::vec3::Vec3;

pub struct TrackPoint {
    pub r: Vec3,         // position
    pub t: f64,          // time
    pub E: Option<f64>,  // energy deposit (`None` if no energy deposited)
}

pub struct Track {
    pub points: Vec<TrackPoint>,
}

impl Track {
    pub fn new(position: Vec3, time: f64) -> Self {
        Track { points: vec![TrackPoint { r: position, t: time, E: None }] }
    }

    pub fn new_empty() -> Self {
        Track { points: vec![] }
    }

    pub fn record(&mut self, position: Vec3, time: f64, energy_deposit: Option<f64>) {
        self.points.push(TrackPoint { r: position, t: time, E: energy_deposit });
    }
}
