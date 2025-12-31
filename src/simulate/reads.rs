use std::fs::create_dir_all;
use std::time::Duration;

use crate::args::App;
use crate::errors::SimError;
use crate::simulate::get_phred_counts;
use crate::utils::mean_error_from_phred;
use flate2::write::GzEncoder;
use flate2::Compression;
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
) -> Result<Vec<u8>, SimError> {
    let mut phreds: Vec<u8> = Vec::with_capacity(read_len);

    let choices: Vec<usize> = (0..128).collect();

    // Choose the first phred from initial probabilities
    let mut rng = rng();
    let p = choices
        .choose_weighted(&mut rng, |i| initial_counts[*i])
        .map_err(|_| SimError::EmptyPhredCounts)?;

    let mut previous = *p as u8;
    phreds.push(previous);

    for _ in 1..read_len {
        let new_phred_counts = transitional_counts[previous as usize];

        let new_phred = choices
            .choose_weighted(&mut rng, |i| new_phred_counts[*i])
            .copied()
            .map_err(|_| SimError::InvalidMarkovState(previous))?;

        phreds.push(new_phred as u8);
        previous = new_phred as u8;
    }

    Ok(phreds)
}

fn load_fasta(template_fasta: &PathBuf) -> Result<Vec<FastaRecord>, SimError> {
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

pub fn simulate_reads(args: App) -> Result<(), SimError> {
    info!("Generating phred counts from fastq files...");
    let (initial_counts, transitional_counts) = get_phred_counts(args.template_fastq)?;

    info!("Loading template fasta...");
    let fasta_records = load_fasta(&args.template_fasta)?;
    let num_fasta_records = fasta_records.len();

    // Number of fastq files to simulate.
    let samples: Vec<usize> = (1..args.num_samples + 1).collect();

    create_dir_all(&args.outdir)
        .map_err(|_| SimError::OutputDirectoryError(args.outdir.clone()))?;

    // Simulate reads in parallel.
    info!("Simulating reads...");
    let spinner = ProgressBar::new_spinner();
    spinner.enable_steady_tick(Duration::from_millis(200));
    spinner.set_style(ProgressStyle::with_template(
        "{spinner:.blue} [{elapsed_precise}]",
    )?);

    samples.par_iter().try_for_each(|sample_index| {
        // Thread specific rng.
        let mut rng = rng();

        // Define output fastq file.
        let mut fastq_out = args.outdir.clone();
        fastq_out.push(format!("sample_{}.fastq.gz", sample_index));

        // We'll write to a different file per thread so this is safe.
        let mut fastq_writer = BufWriter::new(GzEncoder::new(
            File::create(&fastq_out)?,
            Compression::fast(),
        ));

        // Define output read tsv.
        let mut read_tsv = args.outdir.clone();
        read_tsv.push(format!("sample_{}.tsv", sample_index));

        let mut read_tsv_writer = BufWriter::new(File::create(&read_tsv)?);
        read_tsv_writer.write_all(b"read_name\tmean_error\tnum_substitutions\n")?;

        // Randomly choose how many sequences to simulate reads from.
        let num_sequences = args
            .num_sequences
            .choose(&mut rng)
            .copied()
            .unwrap_or(1)
            .min(num_fasta_records);

        // Here, we want to randomize fasta records from provided fasta file.
        let random_fasta_sequences = fasta_records.choose_multiple(&mut rng, num_sequences);

        // For each sequence, we randomly choose how many reads to simulate.
        // Here, we actually want to loop over each set of fasta sequences.
        for fasta_record in random_fasta_sequences {
            let num_reads = args.reads_per_sequence.choose(&mut rng).copied().unwrap_or(1);

            for nread in 1..num_reads + 1 {
                let read_len: usize = fasta_record.sequence.len();
                let read_id = format!("{}|{}", fasta_record.id, nread);

                let phred = simulate_phred(read_len, &initial_counts, &transitional_counts)?;

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

                    // Safe: choices always has 3 elements
                    let new_base = choices.choose(&mut rng).expect("choices array is non-empty");
                    new_sequence[*base_index] = *new_base
                }

                // Write simulated read to fastq file.
                write!(
                    fastq_writer,
                    "@{}\n{}\n+\n{}\n",
                    read_id,
                    String::from_utf8_lossy(&new_sequence),
                    String::from_utf8_lossy(&phred)
                )?;

                // Write read info
                writeln!(
                    read_tsv_writer,
                    "{}\t{}\t{}",
                    read_id, mean_error, expected_num_erroneous_bases
                )?;
            }
        }

        fastq_writer.flush()?;
        read_tsv_writer.flush()?;
        Ok::<(), SimError>(())
    })?;

    spinner.finish();
    Ok(())
}
