use rand::Rng;
use rand_distr::{Normal, Distribution};

use std::ops::{Add, AddAssign, Sub, Mul, Div, Neg};
use std::cmp::{PartialEq};
use approx::relative_eq;
use crate::utils::operations::orthonormal_basis;

#[derive(Copy, Clone, Debug)]
pub struct Vec3(pub f64, pub f64, pub f64);

impl Vec3 {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3(x, y, z)
    }

    pub fn mag(self) -> f64 {
        (self.0*self.0 + self.1*self.1 + self.2*self.2).sqrt()
    }

    pub fn norm(self) -> Vec3 {
        let magnitude = self.mag();
        match magnitude {
            0.0 => Vec3(0.0, 0.0, 0.0),
            _   => self / magnitude,
        }
    }

    pub fn dot(self, rhs: Vec3) -> f64 {
        self.0*rhs.0 + self.1*rhs.1 + self.2*rhs.2
    }

    pub fn cross(self, rhs: Vec3) -> Vec3 {
        Vec3(self.1*rhs.2 - self.2*rhs.1,
             self.2*rhs.0 - self.0*rhs.2,
             self.0*rhs.1 - self.1*rhs.0)
    }

    // Deflect vector by an angle sampled from Gaussian with sigma = theta0
    pub fn deflect(&mut self, rng: &mut impl Rng, theta0: f64) {
        let (u, v) = orthonormal_basis(*self);
        let sigma = theta0 / std::f64::consts::SQRT_2;
        let gaussian = Normal::new(0.0, sigma).unwrap();

        let thetax = gaussian.sample(rng);
        let thetay = gaussian.sample(rng);

        let magnitude = self.mag();
        *self = (*self + thetax * u + thetay * v).norm() * magnitude;
    }

    // Rotate vector by exact angle
    pub fn rotate(&mut self, rng: &mut impl Rng, theta: f64) {
        let (u, v) = orthonormal_basis(*self);
        let phi = rng.random::<f64>() * 2.0 * std::f64::consts::PI;
        let magnitude = self.mag();
        *self = (*self * f64::cos(theta) + u * f64::sin(theta) * f64::cos(phi) + v * f64::sin(theta) * f64::sin(phi)).norm() * magnitude;
    }
}

// Addition (+)
impl Add<Vec3> for Vec3 {
    type Output = Vec3;

    fn add(self, rhs: Vec3) -> Vec3 {
        Vec3(self.0+rhs.0, self.1+rhs.1, self.2+rhs.2)
    }
}

impl AddAssign<Vec3> for Vec3 {
    fn add_assign(&mut self, rhs: Vec3) {
        self.0 = self.0 + rhs.0;
        self.1 = self.1 + rhs.1;
        self.2 = self.2 + rhs.2;
    }
}

// Subtraction (-)
impl Sub<Vec3> for Vec3 {
    type Output = Vec3;

    fn sub(self, rhs: Vec3) -> Vec3 {
        Vec3(self.0-rhs.0, self.1-rhs.1, self.2-rhs.2)
    }
}

// Multiplication (*)
impl Mul<f64> for Vec3 {
    type Output = Vec3;

    fn mul(self, rhs: f64) -> Vec3 {
        Vec3(self.0*rhs, self.1*rhs, self.2*rhs)
    }
}

impl Mul<Vec3> for f64 {
    type Output = Vec3;

    fn mul(self, rhs: Vec3) -> Vec3 {
        Vec3(self*rhs.0, self*rhs.1, self*rhs.2)
    }
}

// Division (/)
impl Div<f64> for Vec3 {
    type Output = Vec3;

    fn div(self, rhs: f64) -> Vec3 {
        Vec3(self.0/rhs, self.1/rhs, self.2/rhs)
    }
}

// Negation (-)
impl Neg for Vec3 {
    type Output = Vec3;

    fn neg(self) -> Vec3 {
        Vec3(-self.0, -self.1, -self.2)
    }
}

// Comparison (==)
impl PartialEq<Vec3> for Vec3 {
    fn eq(&self, rhs: &Vec3) -> bool {
        relative_eq!(self.0, rhs.0) && relative_eq!(self.1, rhs.1) && relative_eq!(self.2, rhs.2)
    }
}


// Tests
#[cfg(test)]
mod tests {
    use super::*;
    use crate::assert_vec3_eq;
    use approx::assert_relative_eq;

    #[test]
    fn test_vec3_indexing() {
        let v1 = Vec3(1.0, 2.0, 3.0);
        let v2 = Vec3(5.5, 3.2, 4.5);
        assert_relative_eq!(v1.0, 1.0);
        assert_relative_eq!(v1.1, 2.0);
        assert_relative_eq!(v1.2, 3.0);
        assert_relative_eq!(v2.0, 5.5);
        assert_relative_eq!(v2.1, 3.2);
        assert_relative_eq!(v2.2, 4.5);
    }

    #[test]
    fn test_vec3_add() {
        let v1 = Vec3(0.0, 0.0, 0.0);
        let v2 = Vec3(1.0, 2.0, 3.0);
        let v3 = Vec3(3.7, -0.8, -1.5);
        assert_vec3_eq!(v1 + v2, Vec3(1.0, 2.0, 3.0));
        assert_vec3_eq!(v1 + v3, Vec3(3.7, -0.8, -1.5));
        assert_vec3_eq!(v2 + v3, Vec3(4.7, 1.2, 1.5));
        assert_vec3_eq!(v1 + v2, v2 + v1);
        assert_vec3_eq!(v1 + v3, v3 + v1);
        assert_vec3_eq!(v2 + v3, v3 + v2);
    }

    #[test]
    fn test_vec3_add_assign() {
        let mut v1 = Vec3(0.0, 0.0, 0.0);
        let mut v2 = Vec3(1.0, 2.0, 3.0);
        let v3 = Vec3(3.7, -0.8, -1.5);
        v2 += v1;
        v1 += v2;
        v2 += v3;
        assert_vec3_eq!(v1, Vec3(1.0, 2.0, 3.0));
        assert_vec3_eq!(v2, Vec3(4.7, 1.2, 1.5));
        assert_vec3_eq!(v3, Vec3(3.7, -0.8, -1.5));
    }

    #[test]
    fn test_vec3_sub() {
        let v1 = Vec3(0.0, 0.0, 0.0);
        let v2 = Vec3(1.0, 2.0, 3.0);
        let v3 = Vec3(3.7, -0.8, -1.5);
        assert_vec3_eq!(v2 - v1, Vec3(1.0, 2.0, 3.0));
        assert_vec3_eq!(v1 - v3, Vec3(-3.7, 0.8, 1.5));
        assert_vec3_eq!(v2 - v3, Vec3(-2.7, 2.8, 4.5));
    }

    #[test]
    fn test_vec3_mul() {
        let v1 = Vec3(0.0, 1.0, 2.0);
        let v2 = Vec3(5.2, 3.1, 4.7);
        let v3 = Vec3(-0.9, 7.2, -2.3);
        let f1 = 4.3;
        let f2 = -1.2;
        let f3 = 0.0;
        assert_vec3_eq!(v1 * f1, Vec3(0.0, 4.3, 8.6));
        assert_vec3_eq!(v3 * f2, Vec3(1.08, -8.64, 2.76));
        assert_vec3_eq!(v1 * f1, f1 * v1);
        assert_vec3_eq!(v1 * f2, f2 * v1);
        assert_vec3_eq!(v1 * f3, f3 * v1);
        assert_vec3_eq!(v2 * f1, f1 * v2);
        assert_vec3_eq!(v2 * f2, f2 * v2);
        assert_vec3_eq!(v2 * f3, f3 * v2);
        assert_vec3_eq!(v3 * f1, f1 * v3);
        assert_vec3_eq!(v3 * f2, f2 * v3);
        assert_vec3_eq!(v3 * f3, f3 * v3);
    }

    #[test]
    fn test_vec3_div() {
        let v1 = Vec3(0.0, 1.0, 2.0);
        let v2 = Vec3(1.6, 7.5, 3.2);
        let v3 = Vec3(-4.2, -6.6, 0.5);
        let f1 = 1.0;
        let f2 = -2.0;
        let f3 = 3.2;
        assert_vec3_eq!(v1 / f1, Vec3(0.0, 1.0, 2.0));
        assert_vec3_eq!(v2 / f1, Vec3(1.6, 7.5, 3.2));
        assert_vec3_eq!(v3 / f1, Vec3(-4.2, -6.6, 0.5));
        assert_vec3_eq!(v1 / f2, Vec3(0.0, -0.5, -1.0));
        assert_vec3_eq!(v2 / f2, Vec3(-0.8, -3.75, -1.6));
        assert_vec3_eq!(v3 / f2, Vec3(2.1, 3.3, -0.25));
        assert_vec3_eq!(v1 / f3, Vec3(0.0, 0.3125, 0.625));
        assert_vec3_eq!(v2 / f3, Vec3(0.5, 2.34375, 1.0));
        assert_vec3_eq!(v3 / f3, Vec3(-1.3125, -2.0625, 0.15625));
    }

    #[test]
    fn test_vec3_neg() {
        let v1 = Vec3(0.0, 0.0, 0.0);
        let v2 = Vec3(5.2, 4.6, 1.1);
        let v3 = Vec3(-7.8, 5.6, -1.3);
        assert_vec3_eq!(-v1, v1);
        assert_vec3_eq!(-v2, Vec3(-5.2, -4.6, -1.1));
        assert_vec3_eq!(-v3, Vec3(7.8, -5.6, 1.3));
        assert_vec3_eq!(-(-v1), v1);
        assert_vec3_eq!(-(-v2), v2);
        assert_vec3_eq!(-(-v3), v3);
    }

    #[test]
    fn test_vec3_mag() {
        let v1 = Vec3(0.0, 0.0, 0.0);
        let v2 = Vec3(3.0, 4.0, 0.0);
        let v3 = Vec3(-5.2, 1.5, -2.2);
        assert_relative_eq!(v1.mag(), 0.0);
        assert_relative_eq!(v2.mag(), 5.0);
        assert_relative_eq!(v3.mag(), 5.842088667591412);
        assert_relative_eq!(v1.mag(), (-v1).mag());
        assert_relative_eq!(v2.mag(), (-v2).mag());
        assert_relative_eq!(v3.mag(), (-v3).mag());
    }

    #[test]
    fn test_vec3_norm() {
        let v1 = Vec3(0.0, 0.0, 0.0);
        let v2 = Vec3(3.0, 4.0, 0.0);
        let v3 = Vec3(-5.2, 1.5, -2.2);
        assert_vec3_eq!(v1.norm(), Vec3(0.0, 0.0, 0.0));
        assert_vec3_eq!(v2.norm(), Vec3(0.6, 0.8, 0.0));
        assert_vec3_eq!(v3.norm(), Vec3(-0.8900926185606606, 0.25675748612326743, -0.37657764631412566));
        assert_vec3_eq!(v1.norm(), -(-v1).norm());
        assert_vec3_eq!(v2.norm(), -(-v2).norm());
        assert_vec3_eq!(v3.norm(), -(-v3).norm());
        assert_vec3_eq!(v1, v1.norm() * v1.mag());
        assert_vec3_eq!(v2, v2.norm() * v2.mag());
        assert_vec3_eq!(v3, v3.norm() * v3.mag());
    }

    #[test]
    fn test_vec3_dot_product() {
        let v1 = Vec3(0.0, 1.0, 2.0);
        let v2 = Vec3(5.0, 3.0, 1.0);
        let v3 = Vec3(4.4, -6.5, 1.2);
        assert_relative_eq!(v1.dot(v2), 5.0);
        assert_relative_eq!(v1.dot(v3), -4.1);
        assert_relative_eq!(v2.dot(v3), 3.7);
        assert_relative_eq!(v1.dot(v2), v2.dot(v1));
        assert_relative_eq!(v1.dot(v3), v3.dot(v1));
        assert_relative_eq!(v2.dot(v3), v3.dot(v2));
        assert_relative_eq!(v1.dot(v1), v1.dot(v1));
        assert_relative_eq!(v2.dot(v2), v2.dot(v2));
        assert_relative_eq!(v3.dot(v3), v3.dot(v3));
    }

    #[test]
    fn test_vec3_cross_product() {
        let v1 = Vec3(0.0, 1.0, 2.0);
        let v2 = Vec3(5.3, -1.3, 8.8);
        let v3 = Vec3(-2.1, -5.3, 0.7);
        assert_vec3_eq!(v1.cross(v2), Vec3(11.4, 10.6, -5.3));
        assert_vec3_eq!(v1.cross(v3), Vec3(11.3, -4.2, 2.1));
        assert_vec3_eq!(v2.cross(v3), Vec3(45.73, -22.19, -30.82));
        assert_vec3_eq!(v1.cross(v1), Vec3(0.0, 0.0, 0.0));
        assert_vec3_eq!(v2.cross(v2), Vec3(0.0, 0.0, 0.0));
        assert_vec3_eq!(v3.cross(v3), Vec3(0.0, 0.0, 0.0));
        assert_vec3_eq!(v1.cross(v2), -v2.cross(v1));
        assert_vec3_eq!(v1.cross(v3), -v3.cross(v1));
        assert_vec3_eq!(v2.cross(v3), -v3.cross(v2));
    }
}
