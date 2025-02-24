use core::f64::consts::PI;

// Geosynchronous orbit definitions
// Requirements for an orbit to be called geosynchronous.

/// lower mean motion for geosynchronous satellites
pub const GEOSYNCHRONOUS_LOWER_MEAN_MOTION: f64 = 0.9;
/// upper mean motion for geosynchronous satellites
pub const GEOSYNCHRONOUS_UPPER_MEAN_MOTION: f64 = 1.1;
/// upper eccentricity for geosynchronous satellites
pub const GEOSYNCHRONOUS_ECCENTRICITY_THRESHOLD: f64 = 0.2;
/// upper inclination for geosynchronous satellites
pub const GEOSYNCHRONOUS_INCLINATION_THRESHOLD_DEGREES: f64 = 70.0;

// Mathematical pub constants
// Mathematical convenience pub constants used by sgp4/sdp4 and related routines.

/// Half of three times PI
pub const THREE_PI_HALF: f64 = (3.0 * PI) / 2.0;
/// two thirds
pub const TWO_THIRD: f64 = 6.666_666_666_666_666E-1;
/// Degrees to Radians conversion factor
pub const DEG_TO_RAD: f64 = PI / 180.0;
/// Radians to Degrees conversion factor
pub const RAD_TO_DEG: f64 = 180.0 / PI;

// Time pub constants
// Constants used for time conversions.

/// Number of minutes per day, XMNPDA in spacetrack report #3
pub const MINUTES_PER_DAY: f64 = 1.44E3;
/// Number of seconds per day
pub const SECONDS_PER_DAY: f64 = 8.6400E4;
/// Difference between libpredict's `predict_julian_date_t` and the julian time format used in some of the internal functions
pub const JULIAN_TIME_DIFF: f64 = 2_444_238.5;

// Physical properties
// General physical properties and definitions.

/// J3 Harmonic (WGS '72), XJ3 in spacetrack report #3
pub const J3_HARMONIC_WGS72: f64 = -2.53881E-6;
/// WGS 84 Earth radius km, XKMPER in spacetrack report #3
pub const EARTH_RADIUS_KM_WGS84: f64 = 6.378_137E3;
/// Flattening factor
pub const FLATTENING_FACTOR: f64 = 3.352_810_664_747_48E-3;
/// Earth rotations per siderial day
pub const EARTH_ROTATIONS_PER_SIDERIAL_DAY: f64 = 1.002_737_909_34;
/// Solar radius in kilometers (IAU 76)
pub const SOLAR_RADIUS_KM: f64 = 6.960_00E5;
/// Astronomical unit in kilometers (IAU 76)
pub const ASTRONOMICAL_UNIT_KM: f64 = 1.495_978_706_91E8;
/// Upper elevation threshold for nautical twilight
pub const NAUTICAL_TWILIGHT_SUN_ELEVATION: f64 = -12.0;
/// Speed of light in vacuum
pub const SPEED_OF_LIGHT: f64 = 299_792_458.0;
/// Angular velocity of Earth in radians per seconds
pub const EARTH_ANGULAR_VELOCITY: f64 = 7.292_115E-5;

// \name Iteration pub constants
//  Constants used in iteration functions like predict_max_elevation(),
//  predict_next_aos() and predict_next_los().

pub const FLT_EPSILON: f64 = 1.192_092_90E-07; // decimal constant

/// Threshold used for comparing lower and upper brackets in `find_max_elevation`
pub const MAXELE_TIME_EQUALITY_THRESHOLD: f64 = FLT_EPSILON;
/// Maximum number of iterations in `find_max_elevation`
pub const MAXELE_MAX_NUM_ITERATIONS: f64 = 10_000.0;

// \name General spacetrack report #3 pub constants
//  These pub constants are also used outside of SGP4/SDP4 code.  The
//  pub constants/symbols are defined on page 76 to page 78 in the report.

/// `k_e = sqrt(Newton's universal gravitational * mass of the Earth)`, in units (earth radii per minutes)^3/2
pub const XKE: f64 = 7.436_691_61E-2;
/// Corresponds to `1/2 * J_2 * a_E^2`. `J_2` is the second gravitational zonal harmonic of Earth, `a_E` is the equatorial radius of Earth.
pub const CK2: f64 = 5.413_079E-4;

// \name Specific spacetrack report #3 pub constants
//  These pub constants are only used by SGP4/SDP4 code.  The constants/symbols are
//  defined on page 76 to page 78 in the report.

/// Shorthand for 10^-6.
pub const E6A: f64 = 1.0E-6;
/// Distance units / Earth radii.
pub const AE: f64 = 1.0;
/// Corresponds to `-3/8 * J_4 * a_E^4`, where `J_4` is the fourth gravitational zonal harmonic of Earth.
pub const CK4: f64 = 6.209_887E-7;
/// Parameter for SGP4/SGP8 density function.
pub const S_DENSITY_PARAM: f64 = 1.012_229;
/// Corresponds to `(q_0 - s)^4` in units (earth radii)^4, where `q_0` and s are parameters for the SGP4/SGP8 density function.
pub const QOMS2T: f64 = 1.880_279E-09;

// Constants in deep space subroutines
// Not defined in spacetrack report #3.
//
// The pub constants might originally be defined in Hujsak (1979) and/or Hujsak
// and Hoots (1977), but this is unavailable on the Internet. Reiterated in F.
// R.  Hoots, P. W. Schumacher and R. A. Glober, "A HISTORY OF ANALYTICAL
// ORBIT MODELING IN THE UNITED STATES SPACE SURVEILLANCE SYSTEM", 2004. Page
// numbers below refer to this article.

/// Solar mean motion `(n_s in units radians/minute, p. 29)`
pub const ZNS: f64 = 1.194_59E-5;
/// Solar perturbation coefficient `(C_s in units radians/minute, p. 29)`
pub const C1SS: f64 = 2.986_479_7E-6;
/// Solar eccentricity `(e_s, p. 29)`
pub const ZES: f64 = 1.675E-2;
/// Lunar mean motion `(n_m in units radians/minute, p. 29)`
pub const ZNL: f64 = 1.583_521_8E-4;
/// Lunar perturbation coefficient `(C_m in units radians/minute, p. 29)`
pub const C1L: f64 = 4.796_806_5E-7;
/// Lunar eccentricity `(e_m, p. 29)`
pub const ZEL: f64 = 5.490E-2;
/// Cosine of the solar inclination `(not defined directly in the paper, but corresponds with cos(I_s) with I_s as the solar inclination on p. 29)`
pub const ZCOSIS: f64 = 9.174_486_7E-1;
/// Sine of the solar inclination `(sin(I_s)`, `I_s on p. 29. See comment above)`
pub const ZSINIS: f64 = 3.978_541_6E-1;
/// Corresponds to `sin(omega_s)` `(omega_s defined on p. 29, no description. See comment above)`
pub const ZSINGS: f64 = -9.808_845_8E-1;
/// Corresponds to `cos(omega_s)` `(omega_s defined on p. 29, no description. See comment above)`
pub const ZCOSGS: f64 = 1.945_905E-1;
/// Constants for one-day resonance conditions, satellite-independent for 1-day period satellites `(Initialization of resonance effects of Earth gravity, Q_22, Q_31 and Q_33, p. 31)`
pub const Q22: f64 = 1.789_167_9E-6;
/// See above
pub const Q31: f64 = 2.146_074_8E-6;
/// See above
pub const Q33: f64 = 2.212_301_5E-7;
/// Constants for secular update for resonance effects of Earth gravity `(G_22, G_32, G_44, G_52 and G_54, p. 36)`
pub const G22: f64 = 5.768_639_6;
/// See above
pub const G32: f64 = 9.524_089_8E-1;
/// See above
pub const G44: f64 = 1.801_499_8;
/// See above
pub const G52: f64 = 1.050_833_0;
/// See above
pub const G54: f64 = 4.410_889_8;
/// Constants for 1/2-day resonance conditions, satellite-independent for 1/2-day period satellites (Initialization for resonance effects of Earth gravity, `sqrt(C_ij^2 + S_ij^2)` where ij = 22, 32, 44, 52 and 54, p. 32)
pub const ROOT22: f64 = 1.789_167_9E-6;
/// See above
pub const ROOT32: f64 = 3.739_379_2E-7;
/// See above
pub const ROOT44: f64 = 7.363_695_3E-9;
/// See above
pub const ROOT52: f64 = 1.142_863_9E-7;
/// See above
pub const ROOT54: f64 = 2.176_580_3E-9;
/// The time-derivative of the Greenwich hour angle in radians per minute (\dot{\theta}, used on p. 36. Not directly defined in report, but values and naming are consistent with this)
pub const THDT: f64 = 4.375_269_1E-3;
