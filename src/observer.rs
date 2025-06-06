use crate::consts::{
    DEG_TO_RAD, EARTH_ANGULAR_VELOCITY, EARTH_RADIUS_KM_WGS84, EARTH_ROTATIONS_PER_SIDERIAL_DAY,
    FLATTENING_FACTOR, JULIAN_TIME_DIFF, MAXELE_MAX_NUM_ITERATIONS, MAXELE_TIME_EQUALITY_THRESHOLD,
    MINUTES_PER_DAY, NAUTICAL_TWILIGHT_SUN_ELEVATION, RAD_TO_DEG, SECONDS_PER_DAY,
};
use crate::geodetic::Geodetic;
use crate::julian_date::{predict_to_julian_double, theta_g_jd};
use crate::math::{asin_clamped, fmod2p, vec3_dot, vec3_length, vec3_sub, Vec3};
use crate::orbit::{
    predict_aos_happens, predict_is_geosynchronous, predict_orbit, OrbitPredictionError,
};
use crate::predict::{
    ObserverElements, Pass, Passes, PredictObservation, PredictObserver, PredictPosition,
};
use crate::sun::predict_observe_sun;
use core::f64::consts::PI;

/// Pass stepping direction used for pass stepping function below.
#[non_exhaustive]
#[derive(PartialEq)]
pub enum StepPassDirection {
    PositiveDirection,
    NegativeDirection,
}

/// Refinement mode, AOS or LOS
#[non_exhaustive]
#[derive(Debug, PartialEq)]
pub enum RefineMode {
    AOS,
    LOS,
}

/// Calculates range, azimuth, elevation and relative velocity from the given observer position.
/// and whether orbit is visible if sun elevation is low enough and the orbit is above the horizon, but still in sunlight.
#[must_use]
pub fn predict_observe_orbit(
    observer: &PredictObserver,
    orbit: &PredictPosition,
) -> PredictObservation {
    let mut observation = observer_calculate(observer, orbit);
    observation.visible = false;

    let sun_obs = predict_observe_sun(observer, orbit.time);
    if !orbit.eclipsed
        && (sun_obs.elevation * 180.0 / PI < NAUTICAL_TWILIGHT_SUN_ELEVATION)
        && (observation.elevation * 180.0 / PI > 0.0)
    {
        observation.visible = true;
    }
    observation
}

/// The procedures `observer_calculate` calculates
#[must_use]
#[allow(clippy::similar_names)]
pub fn observer_calculate(
    observer: &PredictObserver,
    orbit: &PredictPosition,
) -> PredictObservation {
    let jultime = predict_to_julian_double(orbit.time) + JULIAN_TIME_DIFF;

    let mut geodetic = Geodetic {
        lat: observer.latitude,
        lon: observer.longitude,
        alt: observer.altitude / 1000.0,
        theta: 0.0,
    };
    let (obs_pos, obs_vel) = calculate_user_posvel(jultime, &mut geodetic);

    let range = vec3_sub(orbit.position, obs_pos);
    let rgvel = vec3_sub(orbit.velocity, obs_vel);

    let range_length = vec3_length(range);
    let range_rate_length = vec3_dot(range, rgvel) / range_length;

    let theta_dot = 2.0 * PI * EARTH_ROTATIONS_PER_SIDERIAL_DAY / SECONDS_PER_DAY;
    let sin_lat = geodetic.lat.sin();
    let cos_lat = geodetic.lat.cos();
    let sin_theta = geodetic.theta.sin();
    let cos_theta = geodetic.theta.cos();

    let top_s = sin_lat * cos_theta * range.0 + sin_lat * sin_theta * range.1 - cos_lat * range.2;
    let top_e = -sin_theta * range.0 + cos_theta * range.1;
    let top_z = cos_lat * cos_theta * range.0 + cos_lat * sin_theta * range.1 + sin_lat * range.2;

    let top_s_dot = sin_lat * (cos_theta * rgvel.0 - sin_theta * range.0 * theta_dot)
        + sin_lat * (sin_theta * rgvel.1 + cos_theta * range.1 * theta_dot)
        - cos_lat * rgvel.2;
    let top_e_dot = -(sin_theta * rgvel.0 + cos_theta * range.0 * theta_dot)
        + (cos_theta * rgvel.1 - sin_theta * range.1 * theta_dot);

    let top_z_dot = cos_lat
        * (cos_theta * (rgvel.0 + range.1 * theta_dot)
            + sin_theta * (rgvel.1 - range.0 * theta_dot))
        + sin_lat * rgvel.2;

    // Azimut
    let y = -top_e / top_s;
    let mut az = (-top_e / top_s).atan();

    if top_s > 0.0 {
        az += PI;
    }
    if az < 0.0 {
        az += 2.0 * PI;
    }

    // Azimut rate
    let y_dot = -(top_e_dot * top_s - top_s_dot * top_e) / (top_s * top_s);
    let az_dot = y_dot / (1.0 + y * y);

    // Elevation
    let x = top_z / range_length;
    let el = asin_clamped(x);

    // Elevation rate
    let x_dot =
        (top_z_dot * range_length - range_rate_length * top_z) / (range_length * range_length);
    let el_dot = x_dot / (1.0 - x * x).sqrt();
    PredictObservation {
        azimuth: az,
        azimuth_rate: az_dot,
        elevation: el,
        elevation_rate: el_dot,
        range: range_length,
        range_rate: range_rate_length,
        range_x: range.0,
        range_y: range.1,
        range_z: range.2,
        time: orbit.time,
        revolutions: orbit.revolutions,
        visible: false,
    }
}

/// # Errors
///
/// `predict_orbit` can fail
pub fn predict_next_aos(
    oe: &ObserverElements,
    start_utc: f64,
) -> Result<PredictObservation, OrbitPredictionError> {
    let mut curr_time = start_utc;
    let mut time_step;

    let mut orbit = predict_orbit(oe.elements, oe.constants, curr_time)?;
    let mut obs = predict_observe_orbit(oe.observer, &orbit);

    //check whether AOS can happen after specified start time
    if predict_aos_happens(oe.elements, oe.observer.latitude)
        && !predict_is_geosynchronous(oe.elements)
        && !orbit.decayed
    {
        // TODO: Time steps have been found in FindAOS/LOS().
        // Might be based on some pre-existing source, root-finding techniques
        // or something. Find them, and improve readability of the code and so that
        // the mathematical stability of the iteration can be checked.
        // Bisection method, Brent's algorithm? Given a coherent root finding algorithm,
        // can rather have one function for iterating the orbit and then let get_next_aos/los
        // specify bounding intervals for the root finding.

        //skip the rest of the pass if the satellite is currently in range, since we want the _next_ AOS.
        if obs.elevation > 0.0 {
            let los = predict_next_los(oe, curr_time)?;
            curr_time = los.time;
            curr_time += 1.0 / (MINUTES_PER_DAY * 1.0) * 20.0; //skip 20 minutes. LOS might still be within the elevation threshold. (rough quickfix from predict)
            orbit = predict_orbit(oe.elements, oe.constants, curr_time)?;
            obs = predict_observe_orbit(oe.observer, &orbit);
        }

        // iteration until the orbit is roughly in range again, before the satellite pass
        while (obs.elevation * 180.0 / PI < -1.0) || (obs.elevation_rate < 0.0) {
            time_step =
                0.00035 * (obs.elevation * 180.0 / PI * ((orbit.altitude / 8400.0) + 0.46) - 2.0);
            curr_time -= time_step;
            orbit = predict_orbit(oe.elements, oe.constants, curr_time)?;
            obs = predict_observe_orbit(oe.observer, &orbit);
        }

        // fine tune the results until the elevation is within a low enough threshold
        while (obs.elevation * 180.0 / PI).abs() > oe.observer.min_elevation {
            time_step = obs.elevation * 180.0 / PI * orbit.altitude.sqrt() / 530_000.0;
            curr_time -= time_step;
            orbit = predict_orbit(oe.elements, oe.constants, curr_time)?;
            obs = predict_observe_orbit(oe.observer, &orbit);
        }
    }
    Ok(obs)
}

/// Rough stepping through a pass. Uses weird time steps from Predict.
///
/// param `observer` Ground station
/// param `orbital_elements` Orbital elements of satellite
/// param `curr_time` Time from which to start stepping
/// param `direction` Either `POSITIVE_DIRECTION` (step from current time to pass end) or `NEGATIVE_DIRECTION` (step from current time to start of pass). In case of the former, the pass will be stepped until either elevation is negative or the derivative of the elevation is negative
/// return Time for when we have stepped out of the pass
/// # Errors
///
/// `predict_orbit` can fail
pub fn step_pass(
    oe: &ObserverElements,
    mut curr_time: f64,
    direction: &StepPassDirection,
) -> Result<(PredictObservation, f64), OrbitPredictionError> {
    loop {
        let orbit = predict_orbit(oe.elements, oe.constants, curr_time)?;
        let obs = predict_observe_orbit(oe.observer, &orbit);

        // weird time stepping from Predict, but which magically works
        let mut time_step = (obs.elevation - 1.0).cos() * orbit.altitude.sqrt() / 25000.0;
        if ((*direction == StepPassDirection::PositiveDirection) && time_step < 0.0)
            || ((*direction == StepPassDirection::NegativeDirection) && time_step > 0.0)
        {
            time_step = -time_step;
        }

        curr_time += time_step;
        if (obs.elevation < oe.observer.min_elevation * DEG_TO_RAD)
            || ((*direction == StepPassDirection::PositiveDirection) && (obs.elevation_rate <= 0.0))
        {
            return Ok((obs, curr_time));
        }
    }
}

/// # Errors
///
/// `predict_orbit` can fail
pub fn refine_obs_elevation(
    oe: &ObserverElements,
    mut curr_time: f64,
    mode: &RefineMode,
) -> Result<(PredictPosition, PredictObservation, f64), OrbitPredictionError> {
    let time_step = 0.001;
    loop {
        let orbit = predict_orbit(oe.elements, oe.constants, curr_time)?;
        let obs = predict_observe_orbit(oe.observer, &orbit);
        curr_time += time_step;
        if *mode == RefineMode::AOS
            && obs.elevation_rate > 0.0
            && (obs.elevation * RAD_TO_DEG > oe.observer.min_elevation)
        {
            return Ok((orbit, obs, curr_time));
        }
        if *mode == RefineMode::LOS
            && obs.elevation_rate < 0.0
            && (obs.elevation * RAD_TO_DEG < oe.observer.min_elevation)
        {
            return Ok((orbit, obs, curr_time));
        }
    }
}

/// # Errors
///
/// `predict_orbit` can fail
pub fn predict_next_los(
    oe: &ObserverElements,
    start_utc: f64,
) -> Result<PredictObservation, OrbitPredictionError> {
    let mut curr_time = start_utc;
    let mut time_step;

    let mut orbit = predict_orbit(oe.elements, oe.constants, curr_time)?;
    let mut obs = predict_observe_orbit(oe.observer, &orbit);

    // check whether AOS/LOS can happen after specified start time
    if predict_aos_happens(oe.elements, oe.observer.latitude)
        && !predict_is_geosynchronous(oe.elements)
        && !orbit.decayed
    {
        // iteration algorithm from Predict, see comments in predict_next_aos().
        // iterate until next satellite pass

        if obs.elevation < 0.0 {
            let aos = predict_next_aos(oe, curr_time)?;
            curr_time = aos.time;
            orbit = predict_orbit(oe.elements, oe.constants, curr_time)?;
            obs = predict_observe_orbit(oe.observer, &orbit);
        }

        // step through the pass
        (_, curr_time) = step_pass(oe, curr_time, &StepPassDirection::PositiveDirection)?;

        // fine tune to elevation threshold
        loop {
            time_step = obs.elevation * 180.0 / PI * orbit.altitude.sqrt() / 502_500.0;
            curr_time += time_step;
            orbit = predict_orbit(oe.elements, oe.constants, curr_time)?;
            obs = predict_observe_orbit(oe.observer, &orbit);
            if (obs.elevation * 180.0 / PI).abs() <= oe.observer.min_elevation {
                break;
            }
        }
    }
    Ok(obs)
}

/// Convenience function for calculation of derivative of elevation at specific time.
///
/// param `observer` Ground station
/// param `orbital_elements` Orbital elements for satellite
/// param `time` Time
/// return Derivative of elevation at input time
/// # Errors
///
/// `predict_orbit` can fail
pub fn elevation_derivative(oe: &ObserverElements, time: f64) -> Result<f64, OrbitPredictionError> {
    let orbit = predict_orbit(oe.elements, oe.constants, time)?;
    let obs = predict_observe_orbit(oe.observer, &orbit);
    Ok(obs.elevation_rate)
}

/// Find maximum elevation bracketed by input lower and upper time.
///
/// param `observer` Ground station
/// param `orbital_elements` Orbital elements of satellite
/// param `lower_time` Lower time bracket
/// param `upper_time` Upper time bracket
/// return Observation of satellite for maximum elevation between lower and upper time brackets
/// # Errors
///
/// Returns `OrbitPredictionError` on failure of various subroutines
pub fn find_max_elevation(
    oe: &ObserverElements,
    lower_time: f64,
    upper_time: f64,
) -> Result<PredictObservation, OrbitPredictionError> {
    let mut max_ele_time_candidate = (upper_time + lower_time) / 2.0;
    let mut iteration = 0.0;
    let mut lower_time = lower_time;
    let mut upper_time = upper_time;
    while ((lower_time - upper_time).abs() > MAXELE_TIME_EQUALITY_THRESHOLD)
        && (iteration < MAXELE_MAX_NUM_ITERATIONS)
    {
        max_ele_time_candidate = (upper_time + lower_time) / 2.0;

        //calculate derivatives for lower, upper and candidate
        let candidate_deriv = elevation_derivative(oe, max_ele_time_candidate)?;
        let lower_deriv = elevation_derivative(oe, lower_time)?;
        let upper_deriv = elevation_derivative(oe, upper_time)?;

        //check whether derivative has changed sign
        if candidate_deriv * lower_deriv < 0.0 {
            upper_time = max_ele_time_candidate;
        } else if candidate_deriv * upper_deriv < 0.0 {
            lower_time = max_ele_time_candidate;
        } else {
            break;
        }
        iteration += 1.0;
    }

    //prepare output
    let orbit = predict_orbit(oe.elements, oe.constants, max_ele_time_candidate)?;
    Ok(predict_observe_orbit(oe.observer, &orbit))
}

/// passes the user's geodetic position and the time of interest and returns the ECI position and velocity of the observer
///
/// The velocity calculation assumes the geodetic position is stationary relative to the earth's surface
///
/// Reference:  The 1992 Astronomical Almanac, page K11.
pub fn calculate_user_posvel(time: f64, geodetic: &mut Geodetic) -> (Vec3, Vec3) {
    geodetic.theta = fmod2p(theta_g_jd(time) + geodetic.lon); /* LMST */

    let c = 1.0
        / (1.0 + FLATTENING_FACTOR * (FLATTENING_FACTOR - 2.0) * geodetic.lat.sin().powf(2.0))
            .sqrt();

    let sq = (1.0 - FLATTENING_FACTOR).powf(2.0) * c;
    let achcp = (EARTH_RADIUS_KM_WGS84 * c + geodetic.alt) * geodetic.lat.cos();
    let obs_pos = Vec3(
        achcp * geodetic.theta.cos(), /* kilometers */
        achcp * geodetic.theta.sin(),
        (EARTH_RADIUS_KM_WGS84 * sq + geodetic.alt) * geodetic.lat.sin(),
    );
    (
        obs_pos,
        Vec3(
            -EARTH_ANGULAR_VELOCITY * obs_pos.1, /* kilometers/second */
            EARTH_ANGULAR_VELOCITY * obs_pos.0,
            0.0,
        ),
    )
}

/// # Errors
///
/// Returns `OrbitPredictionError` on failure of various subroutines
pub fn get_passes(
    oe: &ObserverElements,
    start_utc: f64,
    stop_utc: f64,
) -> Result<Passes, OrbitPredictionError> {
    let mut orbit = predict_orbit(oe.elements, oe.constants, start_utc)?;
    let mut obs = predict_observe_orbit(oe.observer, &orbit);
    let satellite_el = obs.elevation * RAD_TO_DEG;
    let mut currtime = start_utc;
    let mut passes = vec![];

    if satellite_el.abs() >= oe.observer.min_elevation {
        // Already in a pass, find AOS by going backwards in time
        let (_, real_aos) = step_pass(oe, currtime, &StepPassDirection::NegativeDirection)?;
        currtime = real_aos - 1.0;
    }
    'outer: loop {
        let mut pass = Pass {
            aos: None,
            los: None,
            satellite_position_at_aos: None,
            satellite_position_at_los: None,
            max_elevation: None,
        };
        // rough time step until AOS
        loop {
            if currtime >= stop_utc {
                break 'outer;
            }
            orbit = predict_orbit(oe.elements, oe.constants, currtime)?;
            obs = predict_observe_orbit(oe.observer, &orbit);
            let satellite_el = obs.elevation * RAD_TO_DEG;
            if satellite_el >= oe.observer.min_elevation && obs.elevation_rate > 0.0 {
                currtime -= 1.0;
                let (satpos, observation, _) =
                    refine_obs_elevation(oe, currtime, &RefineMode::AOS)?;
                pass.aos = Some(observation);
                pass.satellite_position_at_aos = Some(satpos);
                currtime += 1.0;
                break;
            }
            currtime += 1.0;
        }
        if pass.aos.is_none() {
            return Err(OrbitPredictionError::MissingAOS);
        }
        // now find LOS
        loop {
            orbit = predict_orbit(oe.elements, oe.constants, currtime)?;
            obs = predict_observe_orbit(oe.observer, &orbit);
            let satellite_el = obs.elevation * RAD_TO_DEG;
            if satellite_el <= oe.observer.min_elevation && obs.elevation_rate < 0.0 {
                currtime -= 1.0;
                let (satpos, observation, _) =
                    refine_obs_elevation(oe, currtime, &RefineMode::LOS)?;
                pass.los = Some(observation);
                pass.satellite_position_at_los = Some(satpos);
                currtime += 1.0;
                break;
            }
            currtime += 1.0;
            if currtime >= stop_utc {
                break;
            }
        }
        if pass.aos.is_some() && pass.los.is_some() {
            let maxel_obs = find_max_elevation(
                oe,
                pass.aos.as_ref().expect("already checked").time,
                pass.los.as_ref().expect("already checked").time,
            )?;
            pass.max_elevation = Some(maxel_obs.elevation);
        }
        passes.push(pass);
    }
    Ok(Passes { passes })
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use crate::observer::{get_passes, DEG_TO_RAD};
    use crate::predict::{ObserverElements, PredictObserver};

    // https://github.com/neuromorphicsystems/sgp4/blob/master/src/tle.rs
    #[allow(clippy::float_cmp)]
    fn assert_eq_f64(first: f64, second: f64) {
        if second == 0.0 {
            assert_eq!(first, 0.0);
        } else {
            assert!((first - second).abs() / second < f64::EPSILON);
        }
    }

    #[test]
    fn test_observer() {
        let elements = sgp4::Elements::from_tle(
            Some("RADARSAT-2".to_string()),
            "1 32382U 07061A   25046.11055603  .00000166  00000+0  81462-4 0  9997".as_bytes(),
            "2 32382  98.5808  54.9834 0001175  88.7553 271.3764 14.29984911896451".as_bytes(),
        )
        .unwrap();
        let observer = PredictObserver {
            name: "TEST".to_string(),
            latitude: 50.0 * DEG_TO_RAD,
            longitude: -130.0 * DEG_TO_RAD,
            altitude: 0.0905,
            min_elevation: 0.0,
        };
        let constants = sgp4::Constants::from_elements(&elements).unwrap();
        let oe = ObserverElements {
            elements: &elements,
            constants: &constants,
            observer: &observer,
        };
        #[allow(clippy::cast_precision_loss)]
        let start_epoch = elements.datetime.and_utc().timestamp() as f64;
        let mut passes = get_passes(&oe, start_epoch, start_epoch + 3.0 * 3600.0)
            .unwrap()
            .passes;

        let second = passes.pop().unwrap();
        let first = passes.pop().unwrap();
        assert_eq_f64(first.aos.as_ref().unwrap().time, 1_739_587_541.817_650_3);
        assert_eq_f64(first.los.as_ref().unwrap().time, 1_739_588_409.207_694_5);
        assert_eq_f64(second.aos.as_ref().unwrap().time, 1_739_593_848.217_693_8);
        assert_eq_f64(second.los.as_ref().unwrap().time, 1_739_594_238.896_644_6);
        assert_eq_f64(
            first.aos.as_ref().unwrap().revolutions,
            89_645.818_476_411_63,
        );
        assert_eq_f64(
            second.aos.as_ref().unwrap().revolutions,
            89_646.862_134_186_38,
        );
    }
}
