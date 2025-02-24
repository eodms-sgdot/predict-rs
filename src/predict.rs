use crate::math::Vec3;
use sgp4::{Constants, Elements};

/// Predicted orbital values for satellite at a given time.
#[derive(Debug)]
pub struct PredictPosition {
    /// Timestamp for last call to `orbit_predict`
    pub time: f64,
    /// Whether the orbit has decayed
    pub decayed: bool,
    /// ECI position in km
    pub position: Vec3,
    /// ECI velocity in km/s
    pub velocity: Vec3,
    /// Latitude in radians, northing/easting
    pub latitude: f64,
    /// Longitude in radians, northing/easting
    pub longitude: f64,
    /// Altitude in km
    pub altitude: f64,
    /// Footprint diameter in km
    pub footprint: f64,
    /// Whether satellite is eclipsed by the earth
    pub eclipsed: bool,
    /// Eclipse depth
    pub eclipse_depth: f64,
    /// Number of revolutions
    pub revolutions: f64,
    /*
    /// Orbital phase (mean anomaly)
    pub phase: f64,
    /// The current number of revolutions around Earth

    /// Current inclination (from xinck within sgp4/sdp4)
    pub inclination: f64,
    /// Current right ascension of the ascending node (from xnodek within sgp4/sdp4)
    pub right_ascension: f64,
    /// Current argument of perigee (from omgadf within sgp4/sdp4)
    pub argument_of_perigee: f64,
    */
}

/// Observation point/ground station (QTH)
#[derive(Debug)]
pub struct PredictObserver {
    /// Observatory name
    pub name: String,
    /// Latitude (WGS84, radians)
    pub latitude: f64,
    /// Longitude (WGS84, radians)
    pub longitude: f64,
    /// Altitude (WGS84, meters)
    pub altitude: f64,
    /// Minimum elevation (WGS84, radians)
    pub min_elevation: f64,
}

/// Data relevant for a relative observation of an orbit or similar with respect to an observation point
#[derive(Debug)]
pub struct PredictObservation {
    /// UTC time
    pub time: f64,
    /// Azimuth angle (rad)
    pub azimuth: f64,
    /// Azimuth angle rate (rad/s)
    pub azimuth_rate: f64,
    /// Elevation angle (rad)
    pub elevation: f64,
    /// Elevation angle rate (rad/s)
    pub elevation_rate: f64,
    /// Range (km)
    pub range: f64,
    /// Range x component
    pub range_x: f64,
    /// Range y component
    pub range_y: f64,
    /// Range z component
    pub range_z: f64,
    /// Range velocity (km/s)
    pub range_rate: f64,
    /// Number of revolutions
    pub revolutions: f64,
    /// Visibility status, whether satellite can be seen by optical means.
    /// The satellite is defined to be visible if:
    /// - The satellite is in sunlight
    /// - The satellite is above the horizon
    /// - The sky is dark enough (sun elevation is below a fixed threshold)
    pub visible: bool,
}

/// A pass of a satellite over an observation point
#[derive(Debug)]
pub struct Pass {
    /// acquisiton of signal
    pub aos: Option<PredictObservation>,
    /// satellite position at acquisiton of signal
    pub satellite_position_at_aos: Option<PredictPosition>,
    /// loss of signal
    pub los: Option<PredictObservation>,
    /// satellite position at loss of signal
    pub satellite_position_at_los: Option<PredictPosition>,
    /// max elevation for this pass
    pub max_elevation: Option<f64>,
}

/// Container for a number of passes
#[derive(Debug)]
pub struct Passes {
    pub passes: Vec<Pass>,
}

/// Container for an observer, the orbital elements of a satellite and the constants
pub struct ObserverElements<'a> {
    pub observer: &'a PredictObserver,
    pub elements: &'a Elements,
    pub constants: &'a Constants,
}
