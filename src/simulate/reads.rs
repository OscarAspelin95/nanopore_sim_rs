use std::fs::create_dir_all;
use std::time::Duration;

use crate::args::App;
use crate::simulate::get_phred_counts;
use crate::utils::mean_error_from_phred;
use anyhow::Result;
use flate2::Compression;
use flate2::write::GzEncoder;
use indicatif::{ProgressBar, ProgressStyle};
use log::info;
use needletail::parse_fastx_file;
use rand::{prelude::*, rng};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

struct FastaRecord {
    id: String,
    sequence: Vec<u8>,
}

fn simulate_phred(
    read_len: usize,
    initial_counts: &[usize; 128],
    transitional_counts: &[[usize; 128]; 128],
) -> Vec<u8> {
    let mut phreds: Vec<u8> = Vec::with_capacity(read_len);

    let choices: Vec<usize> = (0..128).collect();

    // Choose the first phred from initial probabilities
    let mut rng = rng();
    let p = choices
        .choose_weighted(&mut rng, |i| initial_counts[*i])
        .unwrap();

    let mut previous = *p as u8;
    phreds.push(previous);

    for _ in 1..read_len {
        let new_phred_counts = transitional_counts[previous as usize];

        let new_phred = choices
            .choose_weighted(&mut rng, |i| new_phred_counts[*i]).copied()
            .unwrap_or_else(|_| panic!("previous: {}, new_choices: {:?}",
                    previous, new_phred_counts));

        phreds.push(new_phred as u8);
        previous = new_phred as u8;
    }

    phreds
}

fn load_fasta<'a>(template_fasta: &PathBuf) -> Result<Vec<FastaRecord>> {
    let mut fasta_reader = parse_fastx_file(template_fasta)?;

    let mut fasta_records: Vec<FastaRecord> = Vec::new();

    while let Some(record) = fasta_reader.next() {
        let record = match record {
            Ok(record) => record,
            Err(_) => continue,
        };

        fasta_records.push(FastaRecord {
            id: String::from_utf8(record.id().to_vec())?,
            sequence: record.seq().to_vec(),
        });
    }

    Ok(fasta_records)
}

pub fn simulate_reads(args: App) -> Result<()> {
    info!("Generating phred counts from fastq files...");
    let (initial_counts, transitional_counts) = get_phred_counts(args.template_fastq).unwrap();

    info!("Loading template fasta...");
    let fasta_records = load_fasta(&args.template_fasta)?;
    let num_fasta_records = fasta_records.len();

    // Number of fastq files to simulate.
    let samples: Vec<usize> = (1..args.num_samples + 1).collect();

    create_dir_all(&args.outdir).unwrap();

    // Simulate reads in parallel.
    info!("Simulating reads...");
    let spinner = ProgressBar::new_spinner();
    spinner.enable_steady_tick(Duration::from_millis(200));
    spinner.set_style(ProgressStyle::with_template(
        "{spinner:.blue} [{elapsed_precise}]",
    )?);

    samples.par_iter().for_each(|sample_index| {
        // Thread specific rng.
        let mut rng = rng();

        // Define output fastq file.
        let mut fastq_out = args.outdir.clone();
        fastq_out.push(format!("sample_{}.fastq.gz", sample_index));

        // We'll write to a different file per thread so this is safe.
        let mut fastq_writer = BufWriter::new(GzEncoder::new(
            File::create(fastq_out).unwrap(),
            Compression::fast(),
        ));

        // Define output read tsv.
        let mut read_tsv = args.outdir.clone();
        read_tsv.push(format!("sample_{}.tsv", sample_index));

        let mut read_tsv_writer = BufWriter::new(File::create(read_tsv).unwrap());
        read_tsv_writer
            .write_all(b"read_name\tmean_error\tnum_substitutions\n")
            .unwrap();

        // Randomly choose how many sequences to simulate reads from.
        let mut num_sequences = args.num_sequences.choose(&mut rng).unwrap();
        num_sequences = num_sequences.min(&num_fasta_records);

        // Here, we want to randomize fasta records from provided fasta file.
        let random_fasta_sequences = fasta_records.choose_multiple(&mut rng, *num_sequences);

        // For each sequence, we randomly choose how many reads to simulate.
        // Here, we actually want to loop over each set of fasta sequences.
        for fasta_record in random_fasta_sequences {
            let num_reads = args.reads_per_sequence.choose(&mut rng).unwrap();

            for nread in 1..*num_reads + 1 {
                let read_len: usize = fasta_record.sequence.len();
                let read_id = format!("{}|{}", fasta_record.id, nread);

                let phred = simulate_phred(read_len, &initial_counts, &transitional_counts);

                let mean_error = mean_error_from_phred(&phred[..]);

                // How many incorrect bases we expect based on simulated phred.
                let expected_num_erroneous_bases =
                    (mean_error * fasta_record.sequence.len() as f64) as usize;

                // Randomly choose what bases to introduce substitutions for.
                let base_index_choices: Vec<usize> = (0..fasta_record.sequence.len()).collect();
                let base_indices_to_substitute =
                    base_index_choices.choose_multiple(&mut rng, expected_num_erroneous_bases);

                // Unfortunately, we need to clone the sequence.
                let mut new_sequence = fasta_record.sequence.clone();

                for base_index in base_indices_to_substitute {
                    let choices = match fasta_record.sequence[*base_index] {
                        b'A' => [b'C', b'G', b'T'],
                        b'C' => [b'A', b'G', b'T'],
                        b'G' => [b'A', b'C', b'T'],
                        b'T' => [b'A', b'C', b'G'],
                        _ => continue,
                    };

                    let new_base = choices.choose(&mut rng).unwrap();
                    new_sequence[*base_index] = *new_base
                }

                // Write simulated read to fastq file.
                fastq_writer.write_all(b"@").unwrap();
                fastq_writer.write_all(read_id.as_bytes()).unwrap();
                fastq_writer.write_all(b"\n").unwrap();
                fastq_writer.write_all(&fasta_record.sequence).unwrap();
                fastq_writer.write_all(b"\n").unwrap();
                fastq_writer.write_all(&phred[..]).unwrap();
                fastq_writer.write_all(b"\n").unwrap();

                // Write read info
                read_tsv_writer.write_all(read_id.as_bytes()).unwrap();
                read_tsv_writer.write_all(b"\t").unwrap();
                read_tsv_writer
                    .write_all(mean_error.to_string().as_bytes())
                    .unwrap();
                read_tsv_writer.write_all(b"\t").unwrap();
                read_tsv_writer
                    .write_all(expected_num_erroneous_bases.to_string().as_bytes())
                    .unwrap();
                read_tsv_writer.write_all(b"\n").unwrap();
            }
        }

        fastq_writer.flush().unwrap();
        read_tsv_writer.flush().unwrap();
    });

    spinner.finish();
    Ok(())
}
