use crate::utils::vec3::Vec3;

// Returns the two missing vectors to create an orthonormal basis with `vec`
pub fn orthonormal_basis(vec: Vec3) -> (Vec3, Vec3) {
    // Create vector non-colinear with `vec`
    let temp = if vec.norm().0 < 0.999 { // check if vec is colinear with (1.0, 0.0, 0.0)
        Vec3(1.0, 0.0, 0.0)
    } else {
        Vec3(0.0, 1.0, 0.0)
    };
    // Create orthonormal basis (vec, u, v)
    let u = vec.cross(temp).norm();
    let v = vec.cross(u).norm();

    (u, v)
}

// Calculate the 8th-degree "log polynomial"
pub fn log_polynomial(p: f64, cs: Vec<f64>) -> f64 {
    cs.into_iter().enumerate().map(|(n, c)| c * f64::ln(p).powf(n as f64)).sum()
}


// Tests
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use crate::assert_vec3_eq;

    #[test]
    fn test_operations_orthonormal_basis() {
        let v1 = Vec3(1.0, 0.0, 0.0);
        let v2 = Vec3(0.0, 1.0, 0.0);
        let v3 = Vec3(0.0, 0.0, 1.0);
        let v4 = Vec3(2.0, 2.0, 2.0);
        let v5 = Vec3(-3.0, -2.0, 6.0);
        assert_vec3_eq!(orthonormal_basis(v1), (Vec3(0.0, 0.0, 1.0), Vec3(0.0, -1.0, 0.0)));
        assert_vec3_eq!(orthonormal_basis(v2), (Vec3(0.0, 0.0, -1.0), Vec3(-1.0, 0.0, 0.0)));
        assert_vec3_eq!(orthonormal_basis(v3), (Vec3(0.0, 1.0, 0.0), Vec3(-1.0, 0.0, 0.0)));
        assert_vec3_eq!(orthonormal_basis(v4), (Vec3(0.0, 1.0/f64::sqrt(2.0), -1.0/f64::sqrt(2.0)), Vec3(-f64::sqrt(2.0/3.0), 1.0/f64::sqrt(6.0), 1.0/f64::sqrt(6.0))));
        assert_vec3_eq!(orthonormal_basis(v5), (Vec3(0.0, 0.9486832980505138, 0.31622776601683794), Vec3(-0.9035079029052513, 0.1355261854357877, -0.40657855630736306)));
    }

    #[test]
    fn test_operations_log_polynomial() {
        let cs1 = vec![1.97185875, -4.90322067e-01, 5.67984147e-01, -3.78515229e-01, 1.96937857e-01, -6.69875048e-02, 1.30714285e-02, -1.31646064e-03, 5.29555090e-05];
        let cs2 = vec![-4.19708069e+05, 1.09601902e+06, -1.24795701e+06, 8.09428640e+05, -3.27064065e+05, 8.42922022e+04, -1.35291466e+04, 1.23625151e+03, -4.92348947e+01];
        let cs3 = vec![1.13754387e+03, -1.13642381e+03, 4.96588219e+02, -1.23563655e+02, 1.91190645e+01, -1.88126582e+00, 1.14850292e-01, -3.97495919e-03, 5.96940644e-05];
        assert_relative_eq!(log_polynomial(40.0, cs1.clone()), 2.117662466747338);
        assert_relative_eq!(log_polynomial(40.0, cs2.clone()), 10.198838490061462);
        assert_relative_eq!(log_polynomial(40.0, cs3.clone()), 10.07478980542952);
        assert_relative_eq!(log_polynomial(100.0, cs1.clone()), 2.186143318736203);
        assert_relative_eq!(log_polynomial(100.0, cs2.clone()), -580.3696781806648);
        assert_relative_eq!(log_polynomial(100.0, cs3.clone()), 3.2174402419224464);
        assert_relative_eq!(log_polynomial(400.0, cs1), 2.3351969588948975);
        assert_relative_eq!(log_polynomial(400.0, cs2), -187889.94823150337);
        assert_relative_eq!(log_polynomial(400.0, cs3), 1.9963390591128984);
    }
}
