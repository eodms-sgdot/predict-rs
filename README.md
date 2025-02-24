# predict-rs

[![Crates.io](https://img.shields.io/crates/v/predict-rs.svg)](https://crates.io/crates/predict-rs)
[![Docs.rs](https://docs.rs/predict-rs/badge.svg)](https://docs.rs/predict-rs)
[![CI](https://github.com/eodms-sgdot/predict-rs/workflows/CI/badge.svg)](https://github.com/eodms-sgdot/predict-rs/actions)

## Getting Started

To get started, just add to Cargo.toml

```toml
[dependencies]
predict-rs = { version = "0.1" }
```

## Usage

```rust
use predict_rs::consts::{DEG_TO_RAD, RAD_TO_DEG};
use predict_rs::observer::get_passes;
use predict_rs::predict::{ObserverElements, PredictObserver};

fn main() {
    let elements = sgp4::Elements::from_tle(
        Some("NOAA 19".to_owned()),
        "1 33591U 09005A   25055.14587056  .00000347  00000+0  20897-3 0  9990".as_bytes(),
        "2 33591  99.0114 119.2046 0012856 303.3262  56.6679 14.13292610827219".as_bytes(),
    )
    .expect("Error creating sgp4::Elements from a TLE");

    let constants =
        sgp4::Constants::from_elements(&elements).expect("Error creating sgp4::Constants from TLE");

    let observer = PredictObserver {
        name: "MY_SITE".to_string(),
        latitude: 30.5 * DEG_TO_RAD,
        longitude: -73.9 * DEG_TO_RAD,
        altitude: 2.5,
        min_elevation: 0.0,
    };
    let oe = ObserverElements {
        elements: &elements,
        constants: &constants,
        observer: &observer,
    };
    let start = chrono::NaiveDateTime::parse_from_str("2025-02-24 14:00:00", "%Y-%m-%d %H:%M:%S")
        .expect("Error parsing timestamps");
    let start_epoch = start.and_utc().timestamp() as f64;
    let Ok(passes) = get_passes(&oe, start_epoch, start_epoch + 2.5 * 3600.0) else {
        panic!("Could not get any passes");
    };
    for pass in passes.passes {
        let aos = pass.aos.expect("Missing AOS");
        let los = pass.los.expect("Missing LOS");
        let start_datetime =
            chrono::DateTime::from_timestamp_micros((aos.time * 1_000_000.0) as i64)
                .expect("Could not convert AOS to timestamp");
        let stop_datetime =
            chrono::DateTime::from_timestamp_micros((los.time * 1_000_000.0) as i64)
                .expect("Could not convert LOS to timestamp");
        let aos_sat = pass
            .satellite_position_at_aos
            .expect("Missing Satellite Position at AOS");
        let los_sat = pass
            .satellite_position_at_los
            .expect("Missing Satellite Position at LOS");
        let aos_sat_lat = aos_sat.latitude * RAD_TO_DEG;
        let los_sat_lat = los_sat.latitude * RAD_TO_DEG;
        let aos_sat_lon = aos_sat.longitude * RAD_TO_DEG;
        let los_sat_lon = los_sat.longitude * RAD_TO_DEG;
        println!("Orbit: {:7.3}, Start: {}, Stop: {} AOS_AZ: {:.3} LOS_AZ: {:7.3}, MAXEL: {:6.3} SATAOSLAT: {:6.3} SATAOSLON: {:7.3} SATLOSLAT: {:6.3} SATLOSLON: {:7.3}",
                 aos.revolutions,
                 start_datetime,
                 stop_datetime,
                 aos.azimuth*RAD_TO_DEG,
                 los.azimuth*RAD_TO_DEG,
                 pass.max_elevation.expect("Missing Max Elevation")*RAD_TO_DEG,
                 aos_sat_lat,aos_sat_lon,los_sat_lat,los_sat_lon,
                 );
    }
}
```

## Copyright and License

* Copyright 1991-2006 John A. Magliacane (KD2BD)
 
* Copyright 2013- Akademisk radioklubb (LA1K)
 
* Copyright 2013-2015 Knut Magnus Kvamtrø (LA3DPA)

* Copyright 2025 His Majesty the King in Right of Canada,
  as represented by the Minister of Natural Resources

Licensed under

 * GNU GPL Version 2.0
   ([LICENSE](LICENSE)

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.

See [CONTRIBUTING.md](CONTRIBUTING.md).
