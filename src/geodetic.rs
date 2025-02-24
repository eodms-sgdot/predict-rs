use core::f64::consts::{FRAC_PI_2, TAU};

use crate::consts::{EARTH_RADIUS_KM_WGS84, FLATTENING_FACTOR, JULIAN_TIME_DIFF};
use crate::julian_date::theta_g_jd;
use crate::math::{fmod2p, Vec3};

#[derive(Debug)]
pub struct Geodetic {
    pub lat: f64,
    pub lon: f64,
    pub alt: f64,
    pub theta: f64,
}

/// Calculate the Geodetic position of a satellite given a unix epoch time
/// and it's ECI position
#[must_use]
pub fn calculate_lat_lon_alt(time: f64, pos: Vec3) -> Geodetic {
    /* Procedure Calculate_LatLonAlt will calculate the geodetic  */
    /* position of an object given its ECI position pos and time. */
    /* It is intended to be used to determine the ground track of */
    /* a satellite.  The calculations  assume the earth to be an  */
    /* oblate spheroid as defined in WGS '72.                     */

    /* Reference:  The 1992 Astronomical Almanac, page K12. */

    //Convert to julian time:
    let newtime = time + JULIAN_TIME_DIFF;
    let theta = pos.1.atan2(pos.0); /* radians */
    let lon = fmod2p(theta - theta_g_jd(newtime)); /* radians */
    let r = (pos.0.powf(2.0) + pos.1.powf(2.0)).sqrt();
    let e2 = FLATTENING_FACTOR * (2.0 - FLATTENING_FACTOR);
    let mut lat = pos.2.atan2(r); /* radians */
    let mut phi;
    let mut c;
    loop {
        phi = lat;
        c = 1.0 / (1.0 - e2 * phi.sin().powf(2.0)).sqrt();
        lat = (pos.2 + EARTH_RADIUS_KM_WGS84 * c * e2 * phi.sin()).atan2(r);
        // geodetic->lat=atan2(pos[2]+EARTH_RADIUS_KM_WGS84*c*e2*sin(phi),r);
        if (lat - phi).abs() < 1E-10 {
            break;
        }
    }

    let alt = r / lat.cos() - EARTH_RADIUS_KM_WGS84 * c; /* kilometers */

    if lat > FRAC_PI_2 {
        lat -= TAU;
    }
    Geodetic {
        lat,
        lon,
        alt,
        theta,
    }
}
