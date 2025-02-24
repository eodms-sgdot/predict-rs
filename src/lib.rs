//! This crate is a port of the libpredict c library [https://github.com/la1k/libpredict](https://github.com/la1k/libpredict) to Rust
//!
//! It provides methods to calculate passes of a satellite over an observer
//!
//! # Examples
//!
//! Examples can be found in the repository [https://github.com/eodms-sgdot/predict-rs/tree/master/examples](https://github.com/eodms-sgdot/predict-rs/tree/master/examples).
//!

/// Constants used by this crate
pub mod consts;
pub mod geodetic;
/// Julian Date conversions
pub mod julian_date;
/// Misc math functions
pub mod math;
pub mod observer;
pub mod orbit;
pub mod predict;
pub mod sun;
