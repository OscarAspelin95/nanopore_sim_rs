mod args;
mod errors;
mod simulate;
mod utils;

use crate::simulate::simulate_reads;
use args::App;
use clap::Parser;
use log::error;
use simple_logger::SimpleLogger;

fn main() {
    if let Err(e) = SimpleLogger::new().init() {
        eprintln!("Failed to initialize logger: {}", e);
        std::process::exit(1);
    }

    let args: App = App::parse();

    if let Err(e) = simulate_reads(args) {
        error!("Error: {}", e);
        std::process::exit(1);
    }
}
