# nanopore_sim_rs
Simulate full length nanopore amplicon reads using a first order Markov chain like algorithm.

## Requirements
- Linux OS (Ubuntu 24.04.2)
- Rust >= 1.88.0

## Installation
Clone the repository or download the source code. Enter the nanopore_sim_rs directory and run:<br>
`cargo build --release`

The generated binary is available in `target/release/nanopore_sim_rs`.

## Usage
Run with:<br>
`nanopore_sim_rs --template_fastq <reads.fastq.gz> --template-fasta <sequences.fasta> --num-sequences <num_sequences> --reads-per-sequence <reads-per-sequence> --outdir <outdir>`<br>


Optional arguments:
<pre>
<b>--num_samples</b> [10]. How many fastq files to simulate.
</pre>

## Notes
- `--template-fastq` - the fastq file to model errors from. If providing multiple files, the aggregated error profile will be used.

- `--template-fasta` - the fasta file to model reads from. Since sequences are chosen randomly, it is a good idea to dereplicate sequences on species level prior to simulating reads to ensure proper species sampling.

- `--num-sequences` - how many fasta sequences to model reads from for a given sample. If providing multiple values `--num-sequences n1 n2 n3...`, a value will be randomized for each sample.

- `--reads-per-sequence` - how many reads to choose for a given sequence and sample. If providing multiple values `--reads-per-sequence r1 r2 r3 ...`, a value will be randomized.

- `--outdir` - outputs one fastq file per sample with an associated .tsv file containing read names, mean error rates and number of introduced substitutions.
