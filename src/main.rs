mod args;
mod simulate;
mod utils;

use crate::simulate::simulate_reads;
use args::App;
use clap::Parser;
use simple_logger::SimpleLogger;

fn main() {
    SimpleLogger::new().init().unwrap();

    let args: App = App::parse();

    let _ = simulate_reads(args);
}
