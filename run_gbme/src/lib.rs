extern crate clap;
extern crate csv;

use clap::{App, Arg};
//use std::collections::HashMap;
use std::error::Error;
use std::process::Command;
use std::{
    env, fs::DirBuilder, path::{Path, PathBuf},
};

// --------------------------------------------------
// type MetaRecord = HashMap<String, String>;

// --------------------------------------------------
type MyResult<T> = Result<T, Box<Error>>;

// --------------------------------------------------
#[derive(Debug)]
pub struct Config {
    matrix: String,
    metadata: Option<String>,
    bin_dir: Option<String>,
    distance: u32,
    euc_dist_percent: f64,
    num_threads: u32,
    num_scans: u32,
    out_dir: PathBuf,
}

// --------------------------------------------------
pub fn run(config: Config) -> MyResult<()> {
    println!("Using input matrix \"{}\"", &config.matrix);

    let out_dir = &config.out_dir;
    if !out_dir.is_dir() {
        DirBuilder::new().recursive(true).create(&out_dir)?;
    }

    let matrix = config.matrix;
    if !Path::new(&matrix).is_file() {
        let msg = format!("-f \"{}\" is not a file", &matrix);
        return Err(From::from(msg));
    }

    if let Some(metadata_file) = config.metadata {
        //make_metadata_dir(metadata_file, out_dir)?;
        println!("Processing metadata \"{}\"", &metadata_file);
        let make_meta = "make_metadata_dir.py";
        let make_meta_path = match &config.bin_dir {
            Some(d) => Path::new(&d).join(make_meta),
            _ => PathBuf::from(make_meta),
        };

        let meta_dir = config.out_dir.join(PathBuf::from("meta"));

        let mut args = vec![
            "-f".to_string(),
            metadata_file.to_string(),
            "-o".to_string(),
            meta_dir.to_path_buf().to_string_lossy().to_string(),
        ];

        let mut dist_args = if config.distance > 0 {
            vec!["-s".to_string(), config.distance.to_owned().to_string()]
        } else {
            vec![]
        };

        args.append(&mut dist_args);

        let mut euc_dist = if config.euc_dist_percent > 0.0 {
            vec!["-e".to_string(), config.euc_dist_percent.to_string()]
        } else {
            vec![]
        };

        args.append(&mut euc_dist);

        let cmd = Command::new(&make_meta_path).args(args).output();

        if let Err(e) = cmd {
            let msg = format!(
                "Failed to run \"{}\": {}",
                make_meta_path.to_string_lossy(),
                e.to_string()
            );
            return Err(From::from(msg));
        };
    }

    let sna = "sna.r";
    let sna_path = match &config.bin_dir {
        Some(d) => Path::new(&d).join(sna),
        _ => PathBuf::from(sna),
    };

    let status = Command::new(&sna_path)
        .arg("-f")
        .arg(&matrix)
        .arg("-o")
        .arg(&out_dir)
        .arg("-n")
        .arg(&config.num_scans.to_string())
        .status()
        .unwrap();

    return if status.success() {
        Ok(())
    } else {
        let msg = format!("Failed to run \"{}\"", sna_path.to_string_lossy());
        Err(From::from(msg))
    };
}

// --------------------------------------------------
pub fn get_args() -> MyResult<Config> {
    let matches = App::new("Mash All vs All")
        .version("0.1.0")
        .author("Ken Youens-Clark <kyclark@email.arizona.edu")
        .about("Run Mash all-vs-all")
        .arg(
            Arg::with_name("matrix")
                .short("f")
                .long("file")
                .value_name("FILE")
                .help("Distance matrix")
                .required(true),
        )
        .arg(
            Arg::with_name("metadata")
                .short("m")
                .long("metadata")
                .value_name("FILE")
                .help("Metadata file"),
        )
        .arg(
            Arg::with_name("out_dir")
                .short("o")
                .long("out_dir")
                .value_name("DIR")
                .help("Output directory"),
        )
        .arg(
            Arg::with_name("euc_dist_percent")
                .short("e")
                .long("euc_dist_percent")
                .value_name("INT")
                .default_value("0.1")
                .help("Euclidean distance percentage"),
        )
        .arg(
            Arg::with_name("sample_distance")
                .short("d")
                .long("distance")
                .value_name("INT")
                .default_value("1000")
                .help("Min. distance to determine \"near\" samples"),
        )
        .arg(
            Arg::with_name("num_scans")
                .short("s")
                .long("scans")
                .value_name("INT")
                .default_value("100000")
                .help("Number of GBME scans"),
        )
        .arg(
            Arg::with_name("num_threads")
                .short("t")
                .long("threads")
                .value_name("INT")
                .default_value("12")
                .help("Number of threads"),
        )
        .arg(
            Arg::with_name("bin_dir")
                .short("b")
                .long("bin_dir")
                .value_name("DIR")
                .help("Location of binaries"),
        )
        .get_matches();

    let metadata = match matches.value_of("metadata") {
        Some(x) => Some(x.to_string()),
        _ => None,
    };

    let out_dir = match matches.value_of("out_dir") {
        Some(x) => PathBuf::from(x),
        _ => {
            let cwd = env::current_dir()?;
            cwd.join(PathBuf::from("gbme-out"))
        }
    };

    let bin_dir = match matches.value_of("bin_dir") {
        Some(x) => Some(x.to_string()),
        _ => None,
    };

    let distance: u32 = match matches.value_of("sample_distance") {
        Some(x) => match x.trim().parse::<u32>() {
            Ok(n) => n,
            _ => 0,
        },
        _ => 0,
    };

    let num_threads: u32 = match matches.value_of("num_threads") {
        Some(x) => match x.trim().parse() {
            Ok(n) if n > 0 && n < 64 => n,
            _ => 0,
        },
        _ => 0,
    };

    let num_scans: u32 = match matches.value_of("num_scans") {
        Some(x) => match x.trim().parse() {
            Ok(n) => n,
            _ => 0,
        },
        _ => 0,
    };

    let euc_dist_percent = match matches.value_of("euc_dist_percent") {
        Some(x) => match x.trim().parse::<f64>() {
            Ok(n) => {
                if n > 1.0 {
                    n / 100.0
                } else {
                    n
                }
            }
            _ => 0.0,
        },
        _ => 0.0,
    };

    let config = Config {
        bin_dir: bin_dir,
        distance: distance,
        euc_dist_percent: euc_dist_percent,
        num_threads: num_threads,
        num_scans: num_scans,
        out_dir: out_dir,
        metadata: metadata,
        matrix: matches.value_of("matrix").unwrap().to_string(),
    };

    Ok(config)
}

// --------------------------------------------------
// fn make_metadata_dir(meta_file: String, out_dir: String) -> MyResult<()> {
//     let f = match File::open(&meta_file) {
//         Ok(file) => file,
//         Err(e) => {
//             let msg = format!("Can't read {}: {}", meta_file, e.to_string());
//             return Err(From::from(msg));
//         }
//     };
//
//     let meta_dir = config.out_dir.join(PathBuf::from("meta"));
//     if !meta_dir.is_dir() {
//         DirBuilder::new().recursive(true).create(&meta_dir)?;
//     }
//
//     let delimiter = match Path::new(&meta_file).extension() {
//         Some(ext) => match ext.to_str() {
//             Some("csv") => b',',
//             _ => b'\t',
//         },
//         _ => b'\t',
//     };
//
//     let mut rdr = csv::ReaderBuilder::new()
//         .delimiter(delimiter)
//         .from_reader(f);
//
//     for result in rdr.deserialize() {
//         let record: MetaRecord = result?;
//         println!("{:?}", record);
//         let name = record.get("sample_name");
//         let alias = record.get("alias");
//
//         match (name, alias) {
//             (Some(name), Some(alias)) => {
//                 aliases.insert(name.to_string(), alias.to_string());
//                 ()
//             }
//             _ => println!("Missing sample_name or alias"),
//         }
//     }
//     Ok(())
// }

// --------------------------------------------------
// fn haversine(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> f64 {
//     // a = sin²(Δφ/2) + cos φ1 ⋅ cos φ2 ⋅ sin²(Δλ/2)
//     // c = 2 ⋅ atan2( √a, √(1−a) )
//     // distance = r ⋅ c
//     // where φ is latitude, λ is longitude, r is earth’s
//     // radius (mean radius = 6,371km);
//     let r = 6371.00;
//     let lat1_rad = lat1.to_radians();
//     let lat2_rad = lat2.to_radians();
//     let long1_rad = long1.to_radians();
//     let long2_rad = long2.to_radians();
//
//     let a = ((lat2_rad - lat1_rad) / 2.00).sin().powf(2.00)
//         + lat1_rad.cos() * lat2_rad.cos() * ((long2_rad - long1_rad) / 2.00).sin().powf(2.00);
//     let c = 2.00 * ((a).sqrt().atan2((1.00 - a).sqrt()));
//     r * c
// }
