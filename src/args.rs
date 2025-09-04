use clap::Parser;
use std::path::PathBuf;

#[derive(Debug, Parser)]
pub struct App {
    #[clap(long, required=true, num_args = 1.., value_delimiter = ' ')]
    pub template_fastq: Vec<PathBuf>,

    #[clap(
        long,
        required = true,
        help = "Fasta file to randomly select sequences from."
    )]
    pub template_fasta: PathBuf,

    #[clap(
        long,
        default_value_t = 10,
        help = "How many samples (fastq files) to simulate."
    )]
    pub num_samples: usize,

    #[clap(
        short,
        long,
        required=true,
        num_args = 1..,
        value_delimiter = ' ',
        help = "Number of sequences from template fasta to randomly pick for each sample."
    )]
    pub num_sequences: Vec<usize>,

    #[clap(
        short,
        long,
        required = true,
        num_args = 1..,
        value_delimiter = ' ',
        help = "Number of reads to randomly simulate per chosen sequence."
    )]
    pub reads_per_sequence: Vec<usize>,

    #[clap(short, long, required = true)]
    pub outdir: PathBuf,
}
