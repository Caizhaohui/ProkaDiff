//! RA consensus and JC accept/reject rules from published breseq methods.
//! Layer 0: in-memory only (no FASTQ / bowtie2).

#![deny(unsafe_code)]

pub mod jc;
pub mod ra;

pub use jc::{accept_junction, is_candidate_junction, JunctionSupport, SubAlignment};
pub use ra::{call_consensus, BaseObs, ConsensusCall, PileupColumn, RaOptions, Strand};
