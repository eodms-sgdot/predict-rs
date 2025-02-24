use core::f64::consts::PI;

#[derive(Clone, Copy, Debug)]
pub struct Vec3(pub f64, pub f64, pub f64);
pub struct Vec4(pub f64, pub f64, pub f64, pub f64);

#[must_use]
pub fn vec3_length(v: Vec3) -> f64 {
    (v.0 * v.0 + v.1 * v.1 + v.2 * v.2).sqrt()
}

#[must_use]
pub fn vec3_dot(v: Vec3, u: Vec3) -> f64 {
    v.0 * u.0 + v.1 * u.1 + v.2 * u.2
}

#[must_use]
pub fn vec3_mul_scalar(v: Vec3, a: f64) -> Vec3 {
    Vec3(v.0 * a, v.1 * a, v.2 * a)
}

#[must_use]
pub fn vec3_sub(v1: Vec3, v2: Vec3) -> Vec3 {
    Vec3(v1.0 - v2.0, v1.1 - v2.1, v1.2 - v2.2)
}

/// Returns mod 2PI of argument
#[must_use]
pub fn fmod2p(x: f64) -> f64 {
    let mut ret_val = x % (2.0 * PI);
    if ret_val < 0.0 {
        ret_val += 2.0 * PI;
    }
    ret_val
}

#[must_use]
pub fn acos_clamped(arg: f64) -> f64 {
    let arg = arg.clamp(-1.0, 1.0);
    arg.acos()
}

#[must_use]
pub fn asin_clamped(arg: f64) -> f64 {
    let arg = arg.clamp(-1.0, 1.0);
    arg.asin()
}

// https://github.com/RustPython/RustPython
#[must_use]
pub fn modf(x: f64) -> (f64, f64) {
    if !x.is_finite() {
        if x.is_infinite() {
            return (0.0_f64.copysign(x), x);
        } else if x.is_nan() {
            return (x, x);
        }
    }

    (x.fract(), x.trunc())
}
