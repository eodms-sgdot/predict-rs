use chrono::{Datelike, Timelike};
use core::f64::consts::{PI, TAU};
use sgp4::{Constants, Elements};

use crate::consts::{
    AE, CK2, EARTH_RADIUS_KM_WGS84, GEOSYNCHRONOUS_ECCENTRICITY_THRESHOLD,
    GEOSYNCHRONOUS_INCLINATION_THRESHOLD_DEGREES, GEOSYNCHRONOUS_LOWER_MEAN_MOTION,
    GEOSYNCHRONOUS_UPPER_MEAN_MOTION, JULIAN_TIME_DIFF, MINUTES_PER_DAY, SOLAR_RADIUS_KM,
    TWO_THIRD, XKE,
};
use crate::geodetic::calculate_lat_lon_alt;
use crate::geodetic::Geodetic;
use crate::julian_date::{daynum, julian_date_of_epoch, predict_to_julian_double};
use crate::math::{
    acos_clamped, asin_clamped, vec3_dot, vec3_length, vec3_mul_scalar, vec3_sub, Vec3, Vec4,
};
use crate::observer::calculate_user_posvel;
use crate::predict::PredictPosition;
use crate::sun::sun_predict;

#[non_exhaustive]
#[derive(Debug)]
pub enum OrbitPredictionError {
    MissingAOS,
    Propagation(sgp4::Error),
}

#[must_use]
pub fn predict_is_geosynchronous(m: &Elements) -> bool {
    (m.mean_motion >= GEOSYNCHRONOUS_LOWER_MEAN_MOTION)
        && (m.mean_motion <= GEOSYNCHRONOUS_UPPER_MEAN_MOTION)
        && ((m.eccentricity).abs() <= GEOSYNCHRONOUS_ECCENTRICITY_THRESHOLD)
        && ((m.inclination).abs() <= GEOSYNCHRONOUS_INCLINATION_THRESHOLD_DEGREES)
}

#[must_use]
pub fn predict_apogee(m: &Elements) -> f64 {
    let sma = 331.25 * ((1440.0 / m.mean_motion).ln() * (2.0 / 3.0)).exp();
    sma * (1.0 + m.eccentricity) - EARTH_RADIUS_KM_WGS84
}

#[must_use]
#[allow(clippy::similar_names)]
pub fn predict_perigee(m: &Elements) -> f64 {
    let xno = m.mean_motion * TAU / MINUTES_PER_DAY;
    let a1 = (XKE / xno).powf(TWO_THIRD);
    let cosio = (m.inclination * PI / 180.0).cos();
    let theta2 = cosio * cosio;
    let x3thm1 = 3.0 * theta2 - 1.0;
    let eosq = m.eccentricity * m.eccentricity;
    let betao2 = 1.0 - eosq;
    let betao = betao2.sqrt();
    let del1 = 1.5 * CK2 * x3thm1 / (a1 * a1 * betao * betao2);
    let ao = a1 * (1.0 - del1 * (0.5 * TWO_THIRD + del1 * (1.0 + 134.0 / 81.0 * del1)));
    let delo = 1.5 * CK2 * x3thm1 / (ao * ao * betao * betao2);
    let aodp = ao / (1.0 - delo);

    (aodp * (1.0 - m.eccentricity) - AE) * EARTH_RADIUS_KM_WGS84
}

/// Returns true if the satellite pointed to by "x" can ever rise above the horizon of the ground station
#[must_use]
pub fn predict_aos_happens(m: &Elements, latitude: f64) -> bool {
    if m.mean_motion == 0.0 {
        false
    } else {
        let mut lin = m.inclination;

        if lin >= 90.0 {
            lin = 180.0 - lin;
        }
        let apogee = predict_apogee(m);

        ((EARTH_RADIUS_KM_WGS84 / (apogee + EARTH_RADIUS_KM_WGS84)).acos() + (lin * PI / 180.0))
            > latitude.abs()
    }
}

/// This is the stuff we need to do repetitively while tracking
/// This is the old `Calc()` function
/// # Errors
///
/// Can return `OrbitPredictionError`
pub fn predict_orbit(
    elements: &Elements,
    constants: &Constants,
    unix_time: f64,
) -> Result<PredictPosition, OrbitPredictionError> {
    #[allow(clippy::cast_precision_loss)]
    let tsince =
        (unix_time - (elements.datetime.and_utc().timestamp_micros() as f64 / 1_000_000.0)) / 60.0;
    let jultime = predict_to_julian_double(unix_time);
    let prediction = constants
        .propagate(sgp4::MinutesSinceEpoch(tsince))
        .map_err(OrbitPredictionError::Propagation)?;

    let position = Vec3(
        prediction.position[0],
        prediction.position[1],
        prediction.position[2],
    );
    let velocity = Vec3(
        prediction.velocity[0],
        prediction.velocity[1],
        prediction.velocity[2],
    );

    let sat_geo = calculate_lat_lon_alt(jultime, position);

    let latitude = sat_geo.lat;
    let longitude = sat_geo.lon;
    let altitude = sat_geo.alt;

    // Calculate solar position
    let solar_vector = sun_predict(jultime);

    // Find eclipse depth and if sat is eclipsed
    let (eclipse_depth, eclipsed) = is_eclipsed(position, solar_vector);

    // Calculate footprint
    let footprint = 2.0
        * EARTH_RADIUS_KM_WGS84
        * (EARTH_RADIUS_KM_WGS84 / (EARTH_RADIUS_KM_WGS84 + altitude)).acos();

    // Calculate current number of revolutions around Earth
    let temp = TAU / MINUTES_PER_DAY / MINUTES_PER_DAY;
    let epoch_day = f64::from(
        elements.datetime.ordinal() * 86400
            + elements.datetime.hour() * 3600
            + elements.datetime.minute() * 60
            + elements.datetime.second(),
    ) / 86400.0;
    let epoch_year = elements.datetime.year() - 2000;
    let epoch = 1000.0 * f64::from(epoch_year) + epoch_day;
    let jul_epoch = julian_date_of_epoch(epoch);
    let age = jultime - jul_epoch + JULIAN_TIME_DIFF;
    let mean_motion = elements.mean_motion * temp * MINUTES_PER_DAY;
    let mean_anomaly = elements.mean_anomaly * PI / 180.0;
    #[allow(clippy::cast_precision_loss)]
    let revolutions = ((mean_motion * MINUTES_PER_DAY / (PI * 2.0) + age * elements.drag_term)
        * age
        + mean_anomaly / (2.0 * PI))
        + elements.revolution_number as f64;

    //calculate whether orbit is decayed
    let decayed = predict_decayed(elements, jultime);
    Ok(PredictPosition {
        position,
        velocity,
        latitude,
        longitude,
        altitude,
        eclipsed,
        footprint,
        decayed,
        eclipse_depth,
        revolutions,
        time: unix_time,
    })
}

#[must_use]
pub fn predict_decayed(elements: &Elements, time: f64) -> bool {
    let satepoch = daynum(1, 0, elements.datetime.year()) + i64::from(elements.datetime.day());
    let mut has_decayed = false;
    #[allow(clippy::cast_precision_loss)]
    if satepoch as f64
        + ((16.666_666 - elements.mean_motion) / (10.0 * (elements.mean_motion_dot).abs()))
        < time
    {
        has_decayed = true;
    }
    has_decayed
}

/// Calculates if a position is eclipsed.
#[must_use]
pub fn is_eclipsed(pos: Vec3, sol: Vec3) -> (f64, bool) {
    // Determine partial eclipse
    let sd_earth = asin_clamped(EARTH_RADIUS_KM_WGS84 / vec3_length(pos));
    let rho = vec3_sub(sol, pos);
    let sd_sun = asin_clamped(SOLAR_RADIUS_KM / vec3_length(rho));
    let earth = vec3_mul_scalar(pos, -1.0);
    let delta = acos_clamped(vec3_dot(sol, earth) / vec3_length(sol) / vec3_length(earth));
    let depth = sd_earth - sd_sun - delta;

    let is_eclipsed = if sd_earth < sd_sun {
        false
    } else {
        depth >= 0.0
    };
    (depth, is_eclipsed)
}

/// The procedures `calculate_obs` calculates the *topocentric* coordinates of the object with ECI position,
/// {pos}, and velocity, {vel}, from location {geodetic} at a given {time}.
///
/// The function `calculate_obs` returns a Vec4 containingthe azimuth,
/// elevation, range, and range rate (in that order) with units of
/// radians, radians, kilometers, and kilometers/second, respectively.
/// The WGS '72 geoid is used and the effect of atmospheric refraction
/// (under standard temperature and pressure) is incorporated into the
/// elevation calculation; the effect of atmospheric refraction on
/// range and range rate has not yet been quantified.
#[must_use]
pub fn calculate_obs(time: f64, pos: Vec3, vel: Vec3, mut geodetic: Geodetic) -> Vec4 {
    let (obs_pos, obs_vel) = calculate_user_posvel(time, &mut geodetic);

    let range = vec3_sub(pos, obs_pos);
    let rgvel = vec3_sub(vel, obs_vel);

    let range_length = vec3_length(range);

    let sin_lat = geodetic.lat.sin();
    let cos_lat = geodetic.lat.cos();
    let sin_theta = geodetic.theta.sin();
    let cos_theta = geodetic.theta.cos();
    let top_s = sin_lat * cos_theta * range.0 + sin_lat * sin_theta * range.1 - cos_lat * range.2;
    let top_e = -sin_theta * range.0 + cos_theta * range.1;
    let top_z = cos_lat * cos_theta * range.0 + cos_lat * sin_theta * range.1 + sin_lat * range.2;
    let mut azim = (-top_e / top_s).atan(); /* Azimuth */

    if top_s > 0.0 {
        azim += PI;
    }

    if azim < 0.0 {
        azim += 2.0 * PI;
    }

    let elevation = asin_clamped(top_z / range_length);

    Vec4(
        azim,
        elevation,
        range_length, /* Range (kilometers)  */
        /* Range Rate (kilometers/second) */
        vec3_dot(range, rgvel) / vec3_length(range),
    )
}
