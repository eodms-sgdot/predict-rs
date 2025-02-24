use chrono::NaiveDateTime;
use core::f64::consts::PI;

use crate::consts::{EARTH_ROTATIONS_PER_SIDERIAL_DAY, SECONDS_PER_DAY};
use crate::math::modf;

#[must_use]
pub fn get_julian_start_day() -> NaiveDateTime {
    NaiveDateTime::parse_from_str("1979-12-31 00:00:00", "%Y-%m-%d %H:%M:%S")
        .expect("Known good date")
}

#[must_use]
#[allow(clippy::cast_precision_loss)] // f64 is good until the year 2112
pub fn predict_to_julian(input_time: f64) -> f64 {
    //get number of seconds since 1979-12-31 00:00:00 UTC, convert to days
    let jul_start_secs = get_julian_start_day().and_utc().timestamp_micros() as f64 / 1_000_000.0;
    let seconds = input_time - jul_start_secs;
    seconds / SECONDS_PER_DAY
}

#[must_use]
pub fn predict_to_julian_double(time: f64) -> f64 {
    predict_to_julian(time) + ((time % 1.0) / SECONDS_PER_DAY)
}

/// The function `julian_date_of_year` calculates the Julian Date
/// of Day 0.0 of {year}. This function is used to calculate the
/// Julian Date of any date by using `julian_date_of_year`, DOY,
/// and fraction of day.
///
/// Astronomical Formulae for Calculators, Jean Meeus,
/// pages 23-25. Calculate Julian Date of 0.0 Jan year
#[must_use]
#[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
pub fn julian_date_of_year(year: f64) -> f64 {
    let year = year as i64 - 1;
    let i = year / 100;
    let a = i;
    let i = a / 4;
    let b = 2 - a + i;

    let mut i = (365.25 * year as f64) as i64;
    i += (30.600_1 * 14.0) as i64;

    i as f64 + 1_720_994.5 + b as f64
}

/// The function `julian_date_of_epoch` returns the Julian Date of
/// an epoch specified in the format used in the NORAD two-line
/// element sets. It has been modified to support dates beyond
/// the year 1999 assuming that two-digit years in the range 00-56
/// correspond to 2000-2056. Until the two-line element set format
/// is changed, it is only valid for dates through 2056 December 31.
///
/// Modified to support Y2K
/// Valid 1957 through 2056
#[must_use]
pub fn julian_date_of_epoch(epoch: f64) -> f64 {
    let (mut day, mut year) = modf(epoch * 1.0E-3);
    day *= 1.0E3;
    if year < 57.0 {
        year += 2000.0;
    } else {
        year += 1900.0;
    }

    julian_date_of_year(year) + day
}

/// calculate `theta_g_jd` from a julian date
///
/// Reference:  The 1992 Astronomical Almanac, page B6.
#[must_use]
pub fn theta_g_jd(jd: f64) -> f64 {
    let (ut, _) = modf(jd + 0.5);
    let jd = jd - ut;
    let tu = (jd - 2_451_545.0) / 36_525.0;
    let mut gmst = 24_110.548_41 + tu * (8_640_184.812_866 + tu * (0.0931_04 - tu * 6.2E-6));
    gmst = (gmst + SECONDS_PER_DAY * EARTH_ROTATIONS_PER_SIDERIAL_DAY * ut) % SECONDS_PER_DAY;
    2.0 * PI * gmst / SECONDS_PER_DAY
}

#[allow(clippy::cast_possible_truncation)]
#[must_use]
pub fn daynum(mut month: u32, day: u32, mut year: i32) -> i64 {
    if month < 3 {
        year -= 1;
        month += 12;
    }
    if year < 57 {
        year += 100;
    }
    let yy = f64::from(year);
    let mm = f64::from(month);
    let mut dn =
        ((365.25 * (yy - 80.0)).floor() - (19.0 + yy / 100.0).floor() + (4.75 + yy / 400.0) - 16.0)
            .floor() as i64;
    dn += i64::from(day) + 30 * i64::from(month) + (0.6 * mm - 0.3).floor() as i64;
    dn
}
