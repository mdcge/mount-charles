#[macro_export]
macro_rules! assert_vec3_eq {
    ($a:expr, $b:expr) => {
        assert!(
            $a == $b,
            "Vec3 not equal: left = {:?}, right = {:?}",
            $a, $b
        );
    };
}
