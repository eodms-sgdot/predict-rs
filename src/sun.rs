use core::f64::consts::PI;

use crate::consts::{ASTRONOMICAL_UNIT_KM, JULIAN_TIME_DIFF, SECONDS_PER_DAY};
use crate::geodetic::Geodetic;
use crate::julian_date::predict_to_julian_double;
use crate::math::Vec3;
use crate::orbit::calculate_obs;
use crate::predict::{PredictObservation, PredictObserver};

/// The function `delta_et` has been added to allow calculations on the position of the sun.  It provides the difference between UT (approximately the same as UTC) and ET (now referred to as TDT). This function is based on a least squares fit of data from 1950 to 1991 and will need to be updated periodically. Values determined using data from 1950-1991 in the 1990 Astronomical Almanac.  See `DELTA_ET.WQ1` for details.
#[must_use]
pub fn delta_et(year: f64) -> f64 {
    26.465 + 0.747_622 * (year - 1950.0) + 1.886_913 * (2.0 * PI * (year - 1975.0) / 33.0).sin()
}

/// Returns angle in radians from argument in degrees.
#[must_use]
pub fn radians(degrees: f64) -> f64 {
    degrees * PI / 180.0
}

/// Returns angle in degrees from argument in radians.
#[must_use]
pub fn degrees(radians: f64) -> f64 {
    radians * 180.0 / PI
}

#[must_use]
#[allow(clippy::many_single_char_names)] // Try not to mess with original codes
pub fn sun_predict(time: f64) -> Vec3 {
    let jul_utc = time + JULIAN_TIME_DIFF;
    let mjd = jul_utc - 2_415_020.0;
    let year = 1900.0 + mjd / 365.25;
    let t = (mjd + delta_et(year) / SECONDS_PER_DAY) / 36525.0;
    let m = radians(
        358.475_83 + (35_999.049_75 * t % 360.0) - (0.000_150 + 0.000_003_3 * t) * t.sqrt() % 360.0,
    );
    let l = radians(279.696_68 + (36_000.768_92 * t % 360.0) + 0.000_302_5 * t.sqrt() % 360.0);
    let e = 0.016_751_04 - (0.000_041_8 + 0.000_000_126 * t) * t;
    let c = radians(
        (1.919_460 - (0.004_789 + 0.000_014 * t) * t) * m.sin()
            + (0.020_094 - 0.000_100 * t) * (2.0 * m).sin()
            + 0.000_293 * (3.0 * m).sin(),
    );
    let o = radians(259.18 - 1934.142 * t % 360.0);
    let lsa = (l + c - radians(0.00569 - 0.00479 * o.sin())) % (2.0 * PI);
    let nu = (m + c) % (2.0 * PI);
    let mut r = 1.000_000_2 * (1.0 - e.sqrt()) / (1.0 + e * nu.cos());
    let eps = radians(
        23.452_294 - (0.013_012_5 + (0.000_001_64 - 0.000_000_503 * t) * t) * t + 0.00256 * o.cos(),
    );
    r *= ASTRONOMICAL_UNIT_KM;

    Vec3(
        r * lsa.cos(),
        r * lsa.sin() * eps.cos(),
        r * lsa.sin() * eps.sin(),
    )
}

#[must_use]
pub fn predict_observe_sun(observer: &PredictObserver, time: f64) -> PredictObservation {
    // Find sun position
    let jultime = predict_to_julian_double(time) + JULIAN_TIME_DIFF;

    let solar_vector = sun_predict(jultime);

    let geodetic = Geodetic {
        lat: observer.latitude,
        lon: observer.longitude,
        alt: observer.altitude / 1000.0,
        theta: 0.0,
    };

    /* Zero vector for initializations */
    let zero_vector = Vec3(0.0, 0.0, 0.0);

    /* Solar observed azimuth and elevation vector  */
    let solar_set = calculate_obs(jultime, solar_vector, zero_vector, geodetic);

    let sun_azi = solar_set.0;
    let sun_ele = solar_set.1;

    let sun_range = 1.0 + ((solar_set.2 - ASTRONOMICAL_UNIT_KM) / ASTRONOMICAL_UNIT_KM);
    let sun_range_rate = 1000.0 * solar_set.3;

    PredictObservation {
        time,
        azimuth: sun_azi,
        elevation: sun_ele,
        range: sun_range,
        range_rate: sun_range_rate,
        azimuth_rate: 0.0,
        elevation_rate: 0.0,
        range_x: 0.0,
        range_y: 0.0,
        range_z: 0.0,
        revolutions: 0.0,
        visible: false,
    }
}
