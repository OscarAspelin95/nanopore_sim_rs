use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum SimError {
    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),

    #[error("FASTQ/FASTA parsing error: {0}")]
    ParseError(#[from] needletail::errors::ParseError),

    #[error("UTF-8 encoding error: {0}")]
    Utf8Error(#[from] std::string::FromUtf8Error),

    #[error("No valid records found in template FASTQ files")]
    EmptyPhredCounts,

    #[error("Invalid Markov chain state: no valid transitions from Phred {0}")]
    InvalidMarkovState(u8),

    #[error("Template style error: {0}")]
    TemplateError(#[from] indicatif::style::TemplateError),

    #[error("Failed to create output directory: {0}")]
    OutputDirectoryError(PathBuf),
}
