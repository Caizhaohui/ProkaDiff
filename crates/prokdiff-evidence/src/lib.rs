//! RA / MC / JC consensus evidence (clean-room vs published breseq methods).
//!
//! Layer 0 tests are in-memory (no FASTQ / bowtie2). Alignment and BAM I/O are
//! used by the engine (`run_sample`) on compute nodes only.

#![deny(unsafe_code)]

pub mod align;
pub mod engine;
pub mod error;
pub mod fasta;
pub mod jc;
pub mod mc;
pub mod pileup;
pub mod ra;

pub use align::FastqInput;
pub use engine::{run_sample, EngineOptions};
pub use error::EvidenceError;
pub use fasta::{read_reference, write_combined_fasta, FastaRecord};
pub use jc::{accept_junction, is_candidate_junction, JunctionSupport, SubAlignment};
pub use ra::{call_consensus, BaseObs, ConsensusCall, PileupColumn, RaOptions, Strand};
