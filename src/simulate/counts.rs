use anyhow::Result;
use needletail::parse_fastx_file;
use std::path::PathBuf;

const PHRED_OFFSET: u8 = 33;

pub fn get_phred_counts(fastqs: Vec<PathBuf>) -> Result<([usize; 128], [[usize; 128]; 128])> {
    // We assume phred scores lie within printable ASCII characters (values >= 33 before phred offset).
    let mut initial_counts = [0_usize; 128];

    // A 128x128 array with phred counts. We keep track of what phreds come
    // directly after a particular phred. E.g., array[63] is itself an array
    // of counts of all phreds that we have seen directly after 63 (which is "?").
    // E.g., array[63][63] = 100 means that directly after a "?", we have seen another
    // "?" 100 times.
    let mut transitional_counts = [[0_usize; 128]; 128];

    for fastq in fastqs {
        let mut reader = parse_fastx_file(&fastq)?;

        while let Some(record) = reader.next() {
            let record = match record {
                Ok(record) => record,
                Err(_) => continue,
            };

            let qual = record.qual().unwrap();

            // --- Initial probabilities
            let q_begin = qual.first().unwrap();
            initial_counts[*q_begin as usize] += 1;

            // --- Transitional probabilities
            let mut previous: Option<u8> = None;

            for q in qual {
                match previous {
                    None => previous = Some(*q),
                    Some(q_previous) => {
                        transitional_counts[q_previous as usize][*q as usize] += 1;
                        previous = Some(*q)
                    }
                }
            }
        }
    }

    Ok((initial_counts, transitional_counts))
}
