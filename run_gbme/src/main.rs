extern crate run_gbme;

use std::process;

fn main() {
    let config = run_gbme::get_args().expect("Could not get arguments");

    match run_gbme::run(config) {
        Err(e) => {
            println!("Error: {}", e);
            process::exit(1);
        }
        Ok(out_dir) => println!("Done, see output in \"{}\"", &out_dir),
    };
}
